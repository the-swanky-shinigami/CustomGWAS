import pandas as pd
import numpy as np
import requests
import logging
import time
import os
import platform
import datetime

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats

# Handle optional dependencies
try:
    import gseapy as gp
    GSEAPY_AVAILABLE = True
except ImportError:
    GSEAPY_AVAILABLE = False

try:
    import mygene
    MYGENE_AVAILABLE = True
except ImportError:
    MYGENE_AVAILABLE = False


def _downsample_points(x, y, max_points=120000, keep_mask=None, random_seed=42):
    """Downsample dense plot points while preserving priority points via keep_mask."""
    x = np.asarray(x)
    y = np.asarray(y)
    n = len(x)
    if n <= max_points:
        return x, y

    if keep_mask is None:
        keep_mask = np.zeros(n, dtype=bool)
    else:
        keep_mask = np.asarray(keep_mask, dtype=bool)

    keep_idx = np.where(keep_mask)[0]
    rest_idx = np.where(~keep_mask)[0]
    budget_rest = max(0, int(max_points) - len(keep_idx))

    if budget_rest <= 0:
        sel_idx = keep_idx[: int(max_points)]
    else:
        if len(rest_idx) > budget_rest:
            rng = np.random.default_rng(int(random_seed))
            sampled_rest = rng.choice(rest_idx, size=budget_rest, replace=False)
            sel_idx = np.concatenate([keep_idx, sampled_rest])
        else:
            sel_idx = np.concatenate([keep_idx, rest_idx])

    sel_idx = np.sort(sel_idx)
    return x[sel_idx], y[sel_idx]


def _downsample_rank_balanced(x, y, max_points=120000, keep_mask=None, keep_fraction_cap=0.5):
    """Downsample while preserving full rank coverage and a bounded set of priority points."""
    x = np.asarray(x)
    y = np.asarray(y)
    n = len(x)
    if n <= max_points:
        return x, y

    if keep_mask is None:
        keep_mask = np.zeros(n, dtype=bool)
    else:
        keep_mask = np.asarray(keep_mask, dtype=bool)

    keep_idx = np.where(keep_mask)[0]
    # Avoid allowing priority points to consume the entire budget.
    keep_cap = max(1, int(max_points * float(keep_fraction_cap)))
    if len(keep_idx) > keep_cap:
        # For QQ-like curves, larger y values are more extreme and should be preferred.
        order = np.argsort(y[keep_idx])[::-1]
        keep_idx = keep_idx[order[:keep_cap]]

    selected = set(keep_idx.tolist())
    need = max(0, int(max_points) - len(selected))
    if need > 0:
        rank_idx = np.linspace(0, n - 1, num=need, dtype=int)
        selected.update(rank_idx.tolist())

    sel_idx = np.array(sorted(selected), dtype=int)
    if len(sel_idx) > max_points:
        sel_idx = sel_idx[np.linspace(0, len(sel_idx) - 1, num=max_points, dtype=int)]
    return x[sel_idx], y[sel_idx]

def get_species_id(species_name):
    """
    Maps common species names to Taxonomy IDs or Ensembl names.
    """
    species_map = {
        'human': 'human',
        'homo sapiens': 'human',
        'mouse': 'mouse',
        'mus musculus': 'mouse',
        'rat': 'rat',
        'rattus norvegicus': 'rat',
        'rice': '39947', # Oryza sativa
        'oryza sativa': '39947',
        'arabidopsis': '3702',
        'arabidopsis thaliana': '3702',
        'maize': '4577',
        'zea mays': '4577'
    }
    return species_map.get(species_name.lower(), species_name)

def annotate_top_hits(top_loci_df, species='human', dist_kb=10):
    """
    Annotates top SNPs with nearest gene information using MyGene.info API.
    
    Parameters:
    -----------
    top_loci_df : pandas.DataFrame
        DataFrame containing 'chrom' and 'pos' columns.
    species : str
        Common name or taxonomy ID of the species.
    dist_kb : int
        Distance window to look for genes (in kb).
        
    Returns:
    --------
    pandas.DataFrame
        Original DataFrame with 'Nearest_Gene', 'Gene_Name', 'Distance' columns added.
    """
    if not MYGENE_AVAILABLE:
        print("Warning: 'mygene' library not installed. Skipping gene annotation.")
        return top_loci_df

    print(f"Annotating {len(top_loci_df)} loci using MyGene.info (Species: {species})...")
    
    mg = mygene.MyGeneInfo()
    species_id = get_species_id(species)
    
    annotated_df = top_loci_df.copy()
    annotated_df['Nearest_Gene_ID'] = None
    annotated_df['Gene_Symbol'] = None
    annotated_df['Gene_Name'] = None
    annotated_df['Distance_to_Gene'] = None
    
    # MyGene.info query_many is efficient, but works on IDs. 
    # For genomic coordinates, we need to query one by one or use a different endpoint.
    # MyGene.info has a 'query' method that supports genomic ranges.
    
    for idx, row in annotated_df.iterrows():
        chrom = row['chrom']
        pos = row['pos']
        
        # Define query range
        q = f"genomic_pos.chr:{chrom} AND genomic_pos.start:{pos-dist_kb*1000} AND genomic_pos.end:{pos+dist_kb*1000}"
        
        try:
            # Query for genes in this range
            res = mg.query(q, species=species_id, fields='symbol,name,genomic_pos', limit=1)
            
            if res and 'hits' in res and len(res['hits']) > 0:
                hit = res['hits'][0]
                annotated_df.at[idx, 'Nearest_Gene_ID'] = hit.get('_id')
                annotated_df.at[idx, 'Gene_Symbol'] = hit.get('symbol')
                annotated_df.at[idx, 'Gene_Name'] = hit.get('name')
                
                # Calculate distance
                if 'genomic_pos' in hit:
                    g_pos = hit['genomic_pos']
                    # Handle if genomic_pos is a list (multiple locations)
                    if isinstance(g_pos, list):
                        g_pos = g_pos[0]
                        
                    start = g_pos.get('start', pos)
                    end = g_pos.get('end', pos)
                    
                    if pos < start:
                        dist = start - pos
                    elif pos > end:
                        dist = pos - end
                    else:
                        dist = 0 # Inside gene
                        
                    annotated_df.at[idx, 'Distance_to_Gene'] = dist
            
            # Be nice to the API
            time.sleep(0.1)
            
        except Exception as e:
            print(f"Warning: Annotation failed for {chrom}:{pos} - {e}")
            continue
            
    return annotated_df

def perform_pathway_enrichment(gene_list, species='human', gene_sets='KEGG_2021_Human'):
    """
    Performs pathway enrichment analysis using Enrichr via GSEApy.
    """
    if not GSEAPY_AVAILABLE:
        print("Warning: 'gseapy' library not installed. Skipping pathway enrichment.")
        return None

    print(f"Performing pathway enrichment for {len(gene_list)} genes...")
    
    # Clean gene list
    clean_genes = [g for g in gene_list if g and isinstance(g, str)]
    
    if len(clean_genes) < 3:
        print("Not enough genes for enrichment analysis (need at least 3).")
        return None
        
    try:
        # Map species to Enrichr library names if needed
        # GSEApy handles this reasonably well, but we might need to adjust 'gene_sets'
        
        enr = gp.enrichr(
            gene_list=clean_genes,
            gene_sets=gene_sets,
            organism=species, # 'human', 'mouse', 'yeast', 'fly', 'fish', 'worm'
            outdir=None, # Don't write to disk
        )
        
        if enr.results.empty:
            print("No significant pathways found.")
            return None
            
        # Filter for significance
        sig_results = enr.results[enr.results['Adjusted P-value'] < 0.05].sort_values('Adjusted P-value')
        
        return sig_results
        
    except Exception as e:
        print(f"Enrichment analysis failed: {e}")
        return None


# ==============================================================================
# PDF Summary Report Generator
# ==============================================================================

def _add_title_page(pdf, config):
    """Adds a professional title page to the PDF report."""
    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis('off')

    # Title block
    ax.text(0.5, 0.82, 'GWAS Summary Report', fontsize=28, fontweight='bold',
            ha='center', va='center', color='#1a1a2e')
    ax.text(0.5, 0.76, config.get('TRAIT_NAME', 'Unknown Trait'), fontsize=20,
            ha='center', va='center', color='#16213e')

    # Horizontal rule
    ax.axhline(y=0.72, xmin=0.15, xmax=0.85, color='#0f3460', linewidth=2)

    # Info block
    info_lines = [
        ('Species', config.get('SPECIES', 'N/A')),
        ('Trait Type', config.get('TRAIT_TYPE', 'continuous')),
        ('Genotype', os.path.basename(config.get('GENO_PATH', '?'))),
        ('Format', config.get('GENO_FORMAT', '?')),
        ('Samples', str(config.get('n_samples', '?'))),
        ('SNPs (post-QC)', str(config.get('n_snps', '?'))),
        ('Date', datetime.datetime.now().strftime('%d %B %Y  %H:%M')),
        ('Seed', str(config.get('RANDOM_SEED', '?'))),
    ]
    y = 0.64
    for label, value in info_lines:
        ax.text(0.30, y, f'{label}:', fontsize=11, ha='right', va='center',
                fontweight='bold', color='#333')
        ax.text(0.33, y, value, fontsize=11, ha='left', va='center', color='#555')
        y -= 0.035

    # Footer
    ax.text(0.5, 0.08, 'Generated by GWAS Pipeline v3', fontsize=9,
            ha='center', va='center', color='#999', style='italic')
    ax.text(0.5, 0.05, f'Python {platform.python_version()} · NumPy {np.__version__} · '
            f'{platform.platform()}', fontsize=7, ha='center', va='center', color='#bbb')

    pdf.savefig(fig, dpi=100)
    plt.close(fig)


def _add_config_page(pdf, config):
    """Adds a configuration summary page."""
    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis('off')
    ax.text(0.5, 0.96, 'Configuration Summary', fontsize=18, fontweight='bold',
            ha='center', va='top', color='#1a1a2e')

    # Organize config into sections
    sections = {
        'Input / Output': ['GENO_PATH', 'GENO_FORMAT', 'PHENO_PATH', 'COVARIATE_PATH',
                           'TRAIT_NAME', 'TRAIT_TYPE', 'SPECIES', 'RUN_LABEL'],
        'Quality Control': ['SAMPLE_CALLRATE_THRESHOLD', 'SNP_CALLRATE_THRESHOLD',
                            'MAF_THRESHOLD', 'RUN_HWE_FILTER', 'HWE_P_THRESHOLD',
                            'IMPUTATION_METHOD', 'HET_RATE_FILTER', 'HET_RATE_N_SD',
                            'DROP_MISSING_PHENOTYPES', 'REMOVE_PHENOTYPE_OUTLIERS',
                            'OUTLIER_METHOD', 'OUTLIER_THRESHOLD'],
        'LD & PCA': ['LD_PRUNE_BEFORE_PCA', 'LD_PRUNE_R2_THRESHOLD',
                      'LD_PRUNE_WINDOW_SIZE', 'NUM_PCS', 'AUTO_SELECT_PCS'],
        'Models': ['RUN_NAIVE_GLM', 'RUN_PCA_GLM', 'RUN_LMM', 'RUN_LOCO_LMM',
                   'RUN_MLM_STEPWISE', 'SIGNIFICANCE_ALPHA'],
        'Imputation & Memory': ['RUN_BEAGLE_IMPUTATION', 'BEAGLE_JAR',
                                 'FORCE_FLOAT32', 'MAX_MEMORY_GB'],
        'Reproducibility': ['RANDOM_SEED', 'LOG_TO_FILE'],
    }

    y = 0.90
    for section, keys in sections.items():
        if y < 0.06:
            break
        ax.text(0.05, y, section, fontsize=11, fontweight='bold', color='#0f3460')
        y -= 0.025
        for k in keys:
            if k in config and y > 0.04:
                v = config[k]
                # Truncate long paths
                vs = str(v)
                if len(vs) > 55:
                    vs = '...' + vs[-52:]
                ax.text(0.08, y, f'{k}:', fontsize=7.5, ha='left', color='#333', fontfamily='monospace')
                ax.text(0.52, y, vs, fontsize=7.5, ha='left', color='#555', fontfamily='monospace')
                y -= 0.019
        y -= 0.012

    pdf.savefig(fig, dpi=100)
    plt.close(fig)


def _add_qc_page(pdf, qc_report, config=None):
    """Adds a cleanly formatted QC summary page."""
    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis('off')
    ax.text(0.5, 0.96, 'Quality Control Report', fontsize=18, fontweight='bold',
            ha='center', va='top', color='#1a1a2e')

    if not qc_report:
        ax.text(0.5, 0.5, 'No QC report data available.', fontsize=12, ha='center', color='#999')
        pdf.savefig(fig, dpi=100)
        plt.close(fig)
        return

    # --- Extract clean scalar counts from the QC report ---
    step = qc_report.get('step_counts', {})

    # Build a nicely structured table of QC metrics
    rows = []
    rows.append(('Initial Samples', str(qc_report.get('samples_initial', '?'))))
    rows.append(('Initial SNPs', str(qc_report.get('snps_initial', '?'))))

    # Samples removed
    missing_ph = qc_report.get('samples_removed_missing_pheno', [])
    n_missing = len(missing_ph) if isinstance(missing_ph, list) else missing_ph
    miss_str = '0 (Imputed)' if config and not config.get('DROP_MISSING_PHENOTYPES', True) else str(n_missing)
    rows.append(('Samples Removed (missing phenotype)', miss_str))

    outlier_list = qc_report.get('samples_removed_outlier', [])
    n_outlier = len(outlier_list) if isinstance(outlier_list, list) else outlier_list
    out_str = '0 (Disabled)' if config and not config.get('REMOVE_PHENOTYPE_OUTLIERS', True) else str(n_outlier)
    rows.append(('Samples Removed (outlier)', out_str))
    if isinstance(outlier_list, list) and outlier_list:
        names = [str(o.get('Sample', o) if isinstance(o, dict) else o) for o in outlier_list[:5]]
        rows.append(('  Outlier IDs', ', '.join(names) + ('...' if len(outlier_list) > 5 else '')))

    sample_qc_list = qc_report.get('samples_removed_sample_qc', [])
    rows.append(('Samples Removed (call-rate / het)', str(len(sample_qc_list)) if isinstance(sample_qc_list, list) else str(sample_qc_list)))

    # SNPs removed
    rows.append(('SNPs Removed (all missing)', str(len(qc_report.get('snps_removed_all_missing', [])))))
    rows.append(('SNPs Removed (call-rate)', str(qc_report.get('snps_removed_callrate', 0))))
    hwe_str = '0 (Disabled)' if config and not config.get('RUN_HWE_FILTER', True) else str(qc_report.get('snps_removed_hwe', 0))
    rows.append(('SNPs Removed (HWE)', hwe_str))
    rows.append(('SNPs Removed (MAF)', str(qc_report.get('snps_removed_maf', 0))))

    inv_chr = qc_report.get('snps_removed_invalid_chrom', [])
    rows.append(('SNPs Removed (invalid chrom)', str(len(inv_chr)) if isinstance(inv_chr, list) else str(inv_chr)))

    cons_chr0 = qc_report.get('snps_consolidated_chr0', [])
    rows.append(('SNPs Consolidated (chr 0)', str(len(cons_chr0)) if isinstance(cons_chr0, list) else str(cons_chr0)))

    # Final counts
    rows.append(('', ''))  # spacer
    rows.append(('Samples After QC', str(step.get('samples_final', '?'))))
    rows.append(('SNPs After QC', str(step.get('snps_final', '?'))))

    # Render as a clean two-column table
    y = 0.87
    ax.text(0.5, y + 0.02, 'Summary', fontsize=13, fontweight='bold', ha='center', color='#0f3460')
    y -= 0.015
    ax.axhline(y=y, xmin=0.10, xmax=0.90, color='#ccc', linewidth=0.5)
    y -= 0.020
    for label, value in rows:
        if y < 0.05:
            break
        if label == '':
            y -= 0.010
            continue
        indent = label.startswith(' ')
        fs = 8.5 if indent else 9.5
        clr = '#777' if indent else '#333'
        ax.text(0.12, y, label, fontsize=fs, ha='left', fontweight='normal' if indent else 'bold', color=clr)
        ax.text(0.72, y, value, fontsize=fs, ha='left', color='#555')
        y -= 0.028

    pdf.savefig(fig, dpi=100)
    plt.close(fig)


def _add_model_summary_page(pdf, summary_df, lambda_values, best_model, bonferroni_threshold):
    """Adds a model comparison summary page with a table."""
    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis('off')
    ax.text(0.5, 0.96, 'Model Comparison', fontsize=18, fontweight='bold',
            ha='center', va='top', color='#1a1a2e')

    if summary_df is None or summary_df.empty:
        ax.text(0.5, 0.5, 'No model results available.', fontsize=12, ha='center', color='#999')
        pdf.savefig(fig, dpi=100)
        plt.close(fig)
        return

    # Model summary table (kept in the upper half to avoid overlap with lambda block).
    cols = summary_df.columns.tolist()
    cell_text = summary_df.values.tolist()

    table = ax.table(
        cellText=cell_text,
        colLabels=cols,
        loc='upper center',
        cellLoc='center',
        colColours=['#e6f2ff'] * len(cols),
        bbox=[0.06, 0.47, 0.88, 0.38],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1.0, 1.2)

    # Lambda values in two columns when needed.
    if lambda_values:
        ax.text(0.5, 0.42, 'Genomic Inflation (λ)', fontsize=13, fontweight='bold',
                ha='center', color='#1a1a2e')

        items = [(k, v) for k, v in lambda_values.items() if v is not None]
        half = (len(items) + 1) // 2
        left = items[:half]
        right = items[half:]

        y_left = 0.39
        y_right = 0.39
        for model_name, lam in left:
            if lam is not None:
                status = '✓ OK' if 0.9 <= lam <= 1.1 else '⚠ CHECK'
                ax.text(0.10, y_left, f'{model_name}:', fontsize=9, ha='left', color='#333')
                ax.text(0.29, y_left, f'λ = {lam:.4f}  {status}', fontsize=9, ha='left',
                        color='#28a745' if 0.9 <= lam <= 1.1 else '#dc3545')
                y_left -= 0.023

        for model_name, lam in right:
            status = '✓ OK' if 0.9 <= lam <= 1.1 else '⚠ CHECK'
            ax.text(0.55, y_right, f'{model_name}:', fontsize=9, ha='left', color='#333')
            ax.text(0.74, y_right, f'λ = {lam:.4f}  {status}', fontsize=9, ha='left',
                    color='#28a745' if 0.9 <= lam <= 1.1 else '#dc3545')
            y_right -= 0.023

    # Best model recommendation
    if best_model:
        ax.text(0.5, 0.10, f'Recommended Model: {best_model}', fontsize=13,
                fontweight='bold', ha='center', color='#0f3460',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='#e6f2ff', edgecolor='#0f3460'))

    ax.text(0.5, 0.05, f'Bonferroni threshold: {bonferroni_threshold:.2e}', fontsize=9,
            ha='center', color='#777')

    pdf.savefig(fig, dpi=100)
    plt.close(fig)


def _plot_manhattan_for_pdf(results_df, p_col, q_col, title, bonf, alpha,
                            suggestive=1e-5, max_points=120000,
                            rasterize_scatter=True, random_seed=42):
    """Renders a Manhattan plot and returns the figure (for embedding in PDF)."""
    try:
        from natsort import natsorted
    except ImportError:
        natsorted = sorted

    df_p = results_df.dropna(subset=[p_col, 'chrom', 'pos']).copy()
    if df_p.empty:
        return None
    df_p['chrom'] = df_p['chrom'].astype(str)
    df_p['pos'] = pd.to_numeric(df_p['pos'], errors='coerce').fillna(0).astype(int)
    df_p['-log10p'] = -np.log10(df_p[p_col].astype(np.float64) + 1e-300)

    sorted_c = natsorted(df_p['chrom'].unique())
    df_p['chrom'] = pd.Categorical(df_p['chrom'], categories=sorted_c, ordered=True)
    df_p = df_p.sort_values(['chrom', 'pos'])

    W, gap = 1000, 50
    offsets = {}
    cur = 0
    for c in sorted_c:
        offsets[c] = cur
        cur += W + gap
    cmax = df_p.groupby('chrom', observed=True)['pos'].max().replace(0, 1)
    df_p['x'] = df_p.apply(lambda r: offsets[r['chrom']] + r['pos'] / (cmax[r['chrom']] + 1e-9) * W, axis=1)
    colors = plt.cm.tab10(np.linspace(0, 1, len(sorted_c)))
    cmap_d = {c: colors[i] for i, c in enumerate(sorted_c)}

    # Preserve highly significant points while reducing background for PDF size.
    keep_mask = (df_p[p_col].astype(float).values <= float(suggestive))
    x_full = df_p['x'].to_numpy()
    y_full = df_p['-log10p'].to_numpy()
    c_full = df_p['chrom'].astype(str).map(cmap_d).to_numpy()
    n_before = len(x_full)
    x_use, y_use = _downsample_points(
        x_full,
        y_full,
        max_points=max_points,
        keep_mask=keep_mask,
        random_seed=random_seed,
    )
    if len(x_use) < n_before:
        # Rebuild matching color vector for selected points.
        sel_mask = np.zeros(n_before, dtype=bool)
        # Recompute selected indices deterministically from x/y selection.
        # Use value-pair lookup with stable fallback by index order.
        # This avoids changing function signatures across the module.
        xy_map = {}
        for i, (xx, yy) in enumerate(zip(x_full, y_full)):
            xy_map.setdefault((xx, yy), []).append(i)
        sel_idx = []
        for xx, yy in zip(x_use, y_use):
            lst = xy_map.get((xx, yy), [])
            if lst:
                sel_idx.append(lst.pop(0))
        if sel_idx:
            sel_mask[np.array(sel_idx, dtype=int)] = True
            c_use = c_full[sel_mask]
        else:
            c_use = c_full[:len(x_use)]
    else:
        c_use = c_full

    fig, ax = plt.subplots(figsize=(16, 6))
    ax.scatter(
        x_use,
        y_use,
        c=c_use,
        s=10,
        alpha=0.7,
        edgecolors='none',
        rasterized=bool(rasterize_scatter),
    )
    ax.axhline(-np.log10(suggestive), color='blue', ls=':', lw=1, label=f'Suggestive ({suggestive:.0e})')
    ax.axhline(-np.log10(bonf), color='red', ls='--', lw=1.5, label=f'Bonferroni ({bonf:.1e})')
    if q_col in df_p.columns:
        fdr_sub = df_p[df_p[q_col] < alpha]
        if not fdr_sub.empty:
            ax.axhline(-np.log10(fdr_sub[p_col].max()), color='green', ls='-.',
                       lw=1.5, label=f'FDR (q<{alpha})')
    ticks = [(offsets[c] + W / 2) for c in sorted_c]
    ax.set_xticks(ticks)
    ax.set_xticklabels(sorted_c, rotation=45, fontsize=8)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('-log10(P)')
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(axis='y', alpha=0.3)
    fig.tight_layout()
    return fig


def _plot_qq_for_pdf(p_values, title, max_points=120000, rasterize_scatter=True, random_seed=42):
    """Renders a QQ plot and returns the figure (for embedding in PDF)."""
    p = pd.Series(p_values).dropna()
    if p.empty:
        return None
    ps = np.sort(p)
    obs = -np.log10(ps)
    exp = -np.log10(np.arange(1, len(ps) + 1) / (len(ps) + 1))
    lam = np.median(stats.chi2.ppf(1 - ps, 1)) / stats.chi2.ppf(0.5, 1)

    if len(exp) > max_points:
        keep_mask = ps <= 1e-5
        exp, obs = _downsample_rank_balanced(
            exp,
            obs,
            max_points=max_points,
            keep_mask=keep_mask,
            keep_fraction_cap=0.5,
        )

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(exp, obs, s=10, alpha=0.6, color='steelblue', rasterized=bool(rasterize_scatter))
    ax.plot([0, exp.max()], [0, exp.max()], 'r--', lw=1)
    ax.set_xlabel('Expected -log10(P)')
    ax.set_ylabel('Observed -log10(P)')
    ax.set_title(f'{title}\n(λ = {lam:.3f})', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig


def _plot_comparative_qq_for_pdf(qq_data, max_points=120000, rasterize_scatter=True, random_seed=42):
    """Renders a comparative QQ plot for all models and returns the figure."""
    if not qq_data:
        return None
    fig, ax = plt.subplots(figsize=(8, 8))
    colors = plt.cm.Set1(np.linspace(0, 1, max(len(qq_data), 3)))
    for i, (name, pvals) in enumerate(qq_data.items()):
        p = pd.Series(pvals).dropna()
        if p.empty:
            continue
        ps = np.sort(p)
        obs = -np.log10(ps)
        exp = -np.log10(np.arange(1, len(ps) + 1) / (len(ps) + 1))
        lam = np.median(stats.chi2.ppf(1 - ps, 1)) / stats.chi2.ppf(0.5, 1)
        if len(exp) > max_points:
            keep_mask = ps <= 1e-5
            exp, obs = _downsample_rank_balanced(
                exp,
                obs,
                max_points=max_points,
                keep_mask=keep_mask,
                keep_fraction_cap=0.5,
            )
        ax.scatter(exp, obs, s=8, alpha=0.5, color=colors[i], label=f'{name} (λ={lam:.3f})', rasterized=bool(rasterize_scatter))
    xmax = ax.get_xlim()[1]
    ax.plot([0, xmax], [0, xmax], 'k--', lw=1, alpha=0.5)
    ax.set_xlabel('Expected -log10(P)')
    ax.set_ylabel('Observed -log10(P)')
    ax.set_title('Comparative QQ Plot – All Models', fontsize=14, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig


def _plot_volcano_for_pdf(results_df, p_col, beta_col, title, p_threshold=1e-5,
                          max_points=120000, rasterize_scatter=True, random_seed=42):
    """Renders a volcano plot and returns the figure."""
    df = results_df.dropna(subset=[p_col, beta_col]).copy()
    if df.empty:
        return None
    df['-log10p'] = -np.log10(df[p_col].astype(float) + 1e-300)
    sig = df['-log10p'] > -np.log10(p_threshold)

    x = df[beta_col].to_numpy()
    y = df['-log10p'].to_numpy()
    keep_mask = sig.to_numpy()
    x_use, y_use = _downsample_points(x, y, max_points=max_points, keep_mask=keep_mask, random_seed=random_seed)

    # Rebuild significance mask for selected points.
    xy_map = {}
    for i, (xx, yy) in enumerate(zip(x, y)):
        xy_map.setdefault((xx, yy), []).append(i)
    sel_idx = []
    for xx, yy in zip(x_use, y_use):
        lst = xy_map.get((xx, yy), [])
        if lst:
            sel_idx.append(lst.pop(0))
    if sel_idx:
        sig_use = keep_mask[np.array(sel_idx, dtype=int)]
    else:
        sig_use = np.zeros(len(x_use), dtype=bool)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(x_use[~sig_use], y_use[~sig_use], s=8, alpha=0.4, color='grey', label='NS', rasterized=bool(rasterize_scatter))
    ax.scatter(x_use[sig_use], y_use[sig_use], s=12, alpha=0.7, color='red', label='Significant', rasterized=bool(rasterize_scatter))
    ax.axhline(-np.log10(p_threshold), color='blue', ls='--', lw=1,
               label=f'Threshold ({p_threshold:.0e})')
    ax.set_xlabel('Effect Size (Beta)')
    ax.set_ylabel('-log10(P)')
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig


def _add_top_loci_page(pdf, top_loci, target_model):
    """Adds a top loci table page."""
    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis('off')
    ax.text(0.5, 0.96, f'Top Associated Loci ({target_model.upper()})', fontsize=18,
            fontweight='bold', ha='center', va='top', color='#1a1a2e')

    if top_loci is None or top_loci.empty:
        ax.text(0.5, 0.5, 'No significant loci found.', fontsize=12, ha='center', color='#999')
        pdf.savefig(fig, dpi=100)
        plt.close(fig)
        return

    # Select display columns
    display_cols = []
    for c in ['SNP_ID', 'chrom', 'pos', f'p_value_{target_model}',
              f'beta_{target_model}', 'Gene_Symbol', 'Dist_Gene']:
        if c in top_loci.columns:
            display_cols.append(c)

    if not display_cols:
        display_cols = top_loci.columns[:6].tolist()

    tl = top_loci[display_cols].head(25).copy()

    # Format scientific notation for p-values
    for c in tl.columns:
        if 'p_value' in c:
            tl[c] = tl[c].apply(lambda x: f'{x:.2e}' if pd.notna(x) else 'NA')
        elif 'beta' in c or 'se_' in c:
            tl[c] = tl[c].apply(lambda x: f'{x:.4f}' if pd.notna(x) else 'NA')

    # Shorten column names for table
    short_names = {c: c.replace(f'_{target_model}', '').replace('p_value', 'P').replace('beta', 'Beta')
                   .replace('Gene_Symbol', 'Gene').replace('Dist_Gene', 'Dist(bp)') for c in display_cols}
    tl.columns = [short_names.get(c, c) for c in tl.columns]

    cell_text = tl.values.tolist()
    col_labels = tl.columns.tolist()

    table = ax.table(cellText=cell_text, colLabels=col_labels, loc='upper center',
                     cellLoc='center', colColours=['#e6f2ff'] * len(col_labels),
                     bbox=[0.02, 0.10, 0.96, 0.82])
    table.auto_set_font_size(False)
    table.set_fontsize(7.5)
    table.scale(1.0, 1.3)

    # Alternate row shading
    for i, row_cells in enumerate(table.get_celld()):
        if isinstance(row_cells, tuple) and row_cells[0] > 0 and row_cells[0] % 2 == 0:
            table[row_cells].set_facecolor('#f8f9fa')

    n_total = len(top_loci)
    if n_total > 25:
        ax.text(0.5, 0.06, f'Showing top 25 of {n_total} loci. Full list in CSV file.',
                fontsize=9, ha='center', color='#777')

    pdf.savefig(fig, dpi=100)
    plt.close(fig)


def _add_enrichment_page(pdf, enrichment_df):
    """Adds a pathway enrichment results page."""
    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis('off')
    ax.text(0.5, 0.96, 'Pathway Enrichment Analysis', fontsize=18, fontweight='bold',
            ha='center', va='top', color='#1a1a2e')

    if enrichment_df is None or enrichment_df.empty:
        ax.text(0.5, 0.5, 'No significant pathways enriched (or enrichment not run).',
                fontsize=12, ha='center', color='#999')
        pdf.savefig(fig, dpi=100)
        plt.close(fig)
        return

    # Show top 15 pathways as horizontal bar plot
    top = enrichment_df.head(15).copy()
    top['-log10(adj.P)'] = -np.log10(top['Adjusted P-value'].astype(float) + 1e-300)
    top = top.sort_values('-log10(adj.P)')

    # Inset axes for bar plot
    ax_bar = fig.add_axes([0.35, 0.15, 0.58, 0.72])
    colors = plt.cm.YlOrRd(np.linspace(0.3, 0.9, len(top)))[::-1]
    bars = ax_bar.barh(range(len(top)), top['-log10(adj.P)'].values, color=colors)
    ax_bar.set_yticks(range(len(top)))

    # Truncate long pathway names
    names = top['Term'].tolist()
    names = [n[:50] + '...' if len(str(n)) > 50 else str(n) for n in names]
    ax_bar.set_yticklabels(names, fontsize=7)
    ax_bar.set_xlabel('-log10(Adjusted P-value)', fontsize=10)
    ax_bar.set_title('Top Enriched Pathways', fontsize=12, fontweight='bold')
    ax_bar.grid(axis='x', alpha=0.3)

    pdf.savefig(fig, dpi=100)
    plt.close(fig)


def _add_cojo_page(pdf, cojo_results):
    """Adds a COJO conditional analysis results page."""
    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis('off')
    ax.text(0.5, 0.96, 'COJO Conditional Analysis', fontsize=18, fontweight='bold',
            ha='center', va='top', color='#1a1a2e')

    if cojo_results is None or cojo_results.empty:
        ax.text(0.5, 0.5, 'No secondary signals found (or COJO not run).',
                fontsize=12, ha='center', color='#999')
        pdf.savefig(fig, dpi=100)
        plt.close(fig)
        return

    ax.text(0.5, 0.90, f'{len(cojo_results)} secondary signals identified',
            fontsize=12, ha='center', color='#0f3460')

    # Table of secondary signals
    cols = [c for c in ['Locus_Lead', 'Secondary_SNP', 'chrom', 'pos', 'p_conditional', 'beta_cond']
            if c in cojo_results.columns]
    if cols:
        show = cojo_results[cols].head(20).copy()
        for c in show.select_dtypes(include='float').columns:
            show[c] = show[c].apply(lambda x: f'{x:.2e}' if abs(x) < 0.01 else f'{x:.4f}')
        try:
            from tabulate import tabulate
            table_str = tabulate(show.values.tolist(), headers=cols, tablefmt='plain', stralign='right')
        except ImportError:
            table_str = show.to_string(index=False)
        y = 0.83
        for line in table_str.split('\n')[:25]:
            ax.text(0.05, y, line, fontsize=7, fontfamily='monospace', color='#333')
            y -= 0.020

    pdf.savefig(fig, dpi=100)
    plt.close(fig)


def _add_prs_page(pdf, prs_summary):
    """Adds a PRS summary results page with bar chart."""
    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis('off')
    ax.text(0.5, 0.96, 'Polygenic Risk Score (Clump + Threshold)', fontsize=18, fontweight='bold',
            ha='center', va='top', color='#1a1a2e')

    if prs_summary is None or prs_summary.empty:
        ax.text(0.5, 0.5, 'PRS not computed (or RUN_PRS=False).',
                fontsize=12, ha='center', color='#999')
        pdf.savefig(fig, dpi=100)
        plt.close(fig)
        return

    # Table: handle both old (R2_score) and new (R2_train, R2_test) column formats
    cols = [c for c in prs_summary.columns if c in ['P_threshold', 'N_SNPs', 'R2_score', 'R2_train', 'R2_test', 'Mean_PRS', 'Std_PRS']]
    if cols:
        try:
            from tabulate import tabulate
            table_str = tabulate(prs_summary[cols].values.tolist(), headers=cols,
                                 tablefmt='plain', stralign='right', floatfmt='.6f')
        except ImportError:
            table_str = prs_summary[cols].to_string(index=False)
        y = 0.88
        for line in table_str.split('\n')[:15]:
            ax.text(0.1, y, line, fontsize=8, fontfamily='monospace', color='#333')
            y -= 0.022

    # Bar chart of R² vs threshold: detect which metric column exists
    metric_col = None
    if 'R2_test' in prs_summary.columns:
        metric_col = 'R2_test'
    elif 'R2_train' in prs_summary.columns:
        metric_col = 'R2_train'
    elif 'R2_score' in prs_summary.columns:
        metric_col = 'R2_score'
    
    if metric_col:
        valid_prs = prs_summary.dropna(subset=[metric_col])
        if not valid_prs.empty and len(valid_prs) > 0:
            ax_bar = fig.add_axes([0.15, 0.10, 0.70, 0.40])
            thresholds = valid_prs['P_threshold'].astype(str).values
            r2_values = valid_prs[metric_col].values
            bars = ax_bar.bar(range(len(r2_values)), r2_values, color='#0f3460', alpha=0.8)
            ax_bar.set_xticks(range(len(thresholds)))
            ax_bar.set_xticklabels(thresholds, rotation=45, fontsize=8)
            ax_bar.set_xlabel('P-value threshold', fontsize=10)
            ax_bar.set_ylabel(f'R² ({metric_col}) (PRS vs Phenotype)', fontsize=10)
            ax_bar.set_title('PRS Predictive Performance', fontsize=12, fontweight='bold')
            ax_bar.grid(axis='y', alpha=0.3)

    pdf.savefig(fig, dpi=100)
    plt.close(fig)


def _add_gwas_catalog_page(pdf, catalog_results):
    """Adds a GWAS Catalog lookup results page."""
    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis('off')
    ax.text(0.5, 0.96, 'GWAS Catalog Cross-Reference', fontsize=18, fontweight='bold',
            ha='center', va='top', color='#1a1a2e')

    if catalog_results is None or catalog_results.empty:
        ax.text(0.5, 0.5, 'No GWAS Catalog matches found (or lookup not run / non-human species).',
                fontsize=12, ha='center', color='#999')
        pdf.savefig(fig, dpi=100)
        plt.close(fig)
        return

    ax.text(0.5, 0.90, f'{len(catalog_results)} catalog associations found near top hits',
            fontsize=12, ha='center', color='#0f3460')

    cols = [c for c in ['Query_SNP', 'Catalog_RSid', 'Catalog_Trait', 'Catalog_P']
            if c in catalog_results.columns]
    if cols:
        show = catalog_results[cols].head(25).copy()
        # Truncate long trait names
        if 'Catalog_Trait' in show.columns:
            show['Catalog_Trait'] = show['Catalog_Trait'].apply(
                lambda x: str(x)[:40] + '...' if len(str(x)) > 40 else str(x))
        try:
            from tabulate import tabulate
            table_str = tabulate(show.values.tolist(), headers=cols, tablefmt='plain', stralign='left')
        except ImportError:
            table_str = show.to_string(index=False)
        y = 0.83
        for line in table_str.split('\n')[:30]:
            ax.text(0.05, y, line, fontsize=7, fontfamily='monospace', color='#333')
            y -= 0.020

    pdf.savefig(fig, dpi=100)
    plt.close(fig)


def _add_session_page(pdf, config, elapsed_seconds=None, peak_memory=None):
    """Adds a session info / footer page."""
    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis('off')
    ax.text(0.5, 0.96, 'Session Information', fontsize=18, fontweight='bold',
            ha='center', va='top', color='#1a1a2e')

    y = 0.86
    info_items = [
        ('Platform', platform.platform()),
        ('Python', platform.python_version()),
        ('NumPy', np.__version__),
        ('Pandas', pd.__version__),
    ]
    try:
        import scipy
        info_items.append(('SciPy', scipy.__version__))
    except:
        pass
    try:
        import statsmodels as _sm
        info_items.append(('Statsmodels', _sm.__version__))
    except:
        pass
    try:
        import plotly
        info_items.append(('Plotly', plotly.__version__))
    except:
        pass

    info_items.append(('Random Seed', str(config.get('RANDOM_SEED', '?'))))
    info_items.append(('Float dtype', 'float32' if config.get('FORCE_FLOAT32') else 'float64'))
    info_items.append(('Max memory', f"{config.get('MAX_MEMORY_GB', '?')} GB"))

    if elapsed_seconds is not None:
        m, s = divmod(int(elapsed_seconds), 60)
        info_items.append(('Runtime', f'{m}m {s}s'))
    if peak_memory is not None:
        info_items.append(('Peak traced memory', peak_memory))

    for label, value in info_items:
        ax.text(0.25, y, f'{label}:', fontsize=10, ha='right', fontweight='bold', color='#333')
        ax.text(0.28, y, str(value), fontsize=10, ha='left', color='#555')
        y -= 0.030

    # Config hash
    cfg_str = f"{config.get('GENO_PATH','')}|{config.get('PHENO_PATH','')}|" \
              f"{config.get('COVARIATE_PATH','')}|{config.get('MAF_THRESHOLD','')}|" \
              f"{config.get('RANDOM_SEED','')}"
    import hashlib
    chash = hashlib.md5(cfg_str.encode()).hexdigest()[:12]
    ax.text(0.25, y, 'Config hash:', fontsize=10, ha='right', fontweight='bold', color='#333')
    ax.text(0.28, y, chash, fontsize=10, ha='left', color='#555', fontfamily='monospace')
    y -= 0.06

    # Output files listing
    output_dir = config.get('OUTPUT_DIR', '')
    if output_dir and os.path.isdir(output_dir):
        ax.text(0.5, y, 'Output Files', fontsize=12, fontweight='bold', ha='center', color='#0f3460')
        y -= 0.030
        for fname in sorted(os.listdir(output_dir)):
            fpath = os.path.join(output_dir, fname)
            sz = os.path.getsize(fpath) if os.path.isfile(fpath) else 0
            if sz < 1024:
                sz_str = f'{sz} B'
            elif sz < 1024 ** 2:
                sz_str = f'{sz / 1024:.1f} KB'
            else:
                sz_str = f'{sz / 1024 ** 2:.1f} MB'
            ax.text(0.15, y, f'  {fname}', fontsize=8, ha='left', color='#555', fontfamily='monospace')
            ax.text(0.85, y, sz_str, fontsize=8, ha='right', color='#999')
            y -= 0.020
            if y < 0.05:
                break

    pdf.savefig(fig, dpi=100)
    plt.close(fig)


def generate_summary_report(report_data, output_path=None, pdf_path=None):
    """
    Generates a comprehensive multi-page PDF report with all GWAS results, plots, and tables.

    Parameters
    ----------
    report_data : dict
        Dictionary containing all pipeline results:
            - config        : dict of all configuration parameters
            - results_df    : full GWAS results DataFrame
            - qc_report     : dict of QC step results
            - summary_df    : model comparison DataFrame
            - lambda_values : dict {model_name: lambda_gc}
            - qq_data       : dict {model_name: p_values_array}
            - top_loci      : DataFrame of top loci
            - enrichment_df : DataFrame from pathway enrichment (or None)
            - best_model    : str name of recommended model
            - target_model  : str short key of target model
            - bonferroni    : float bonferroni threshold
            - model_methods : list of model short keys run
            - elapsed       : float seconds elapsed (optional)
            - peak_memory   : str formatted peak memory (optional)
    output_path : str, optional
        Path for a markdown summary (legacy, still generated).
    pdf_path : str, optional
        Path for the PDF report. If None, derived from output_path.

    Returns
    -------
    str : path to the generated PDF file.
    """
    config       = report_data.get('config', {})
    results_df   = report_data.get('results_df', pd.DataFrame())
    qc_report    = report_data.get('qc_report', {})
    summary_df   = report_data.get('summary_df', pd.DataFrame())
    lambda_values = report_data.get('lambda_values', {})
    qq_data      = report_data.get('qq_data', {})
    top_loci     = report_data.get('top_loci', pd.DataFrame())
    enrichment_df = report_data.get('enrichment_df', None)
    best_model   = report_data.get('best_model', None)
    target_model = report_data.get('target_model', 'lmm')
    bonf         = report_data.get('bonferroni', 5e-8)
    model_methods = report_data.get('model_methods', [])
    alpha        = config.get('SIGNIFICANCE_ALPHA', 0.05)
    report_pdf_dpi = int(config.get('REPORT_PDF_DPI', 72))
    report_pdf_compression = int(config.get('REPORT_PDF_COMPRESSION', 9))
    report_max_points = int(config.get('REPORT_MAX_POINTS', 120000))
    report_rasterize_scatter = bool(config.get('REPORT_RASTERIZE_SCATTER', True))
    report_seed = int(config.get('RANDOM_SEED', 42))
    cojo_results  = report_data.get('cojo_results', pd.DataFrame())
    catalog_results = report_data.get('catalog_results', pd.DataFrame())
    prs_summary   = report_data.get('prs_summary', pd.DataFrame())

    # Determine PDF path
    if pdf_path is None:
        if output_path:
            pdf_path = output_path.replace('.md', '.pdf').replace('.txt', '.pdf')
            if not pdf_path.endswith('.pdf'):
                pdf_path = pdf_path + '.pdf'
        else:
            pdf_path = 'GWAS_Summary_Report.pdf'

    print(f"Generating PDF report -> {pdf_path}")

    # Temporarily switch to Agg backend for PDF rendering, then restore
    orig_backend = matplotlib.get_backend()
    plt.switch_backend('Agg')
    _orig_pdf_comp = matplotlib.rcParams.get('pdf.compression', 6)
    matplotlib.rcParams['pdf.compression'] = report_pdf_compression

    try:
        with PdfPages(pdf_path) as pdf:
            # --- Page 1: Title ---
            _add_title_page(pdf, config)

            # --- Page 2: Configuration ---
            _add_config_page(pdf, config)

            # --- Page 3: QC Summary ---
            _add_qc_page(pdf, qc_report, config)

            # --- Page 4: Model Comparison ---
            _add_model_summary_page(pdf, summary_df, lambda_values, best_model, bonf)

            # --- Pages 5+: Manhattan plots ---
            for method in model_methods:
                p_col = f'p_value_{method}'
                q_col = f'q_value_{method}'
                if p_col in results_df.columns:
                    fig = _plot_manhattan_for_pdf(results_df, p_col, q_col,
                                                  f'Manhattan Plot – {method.upper()}',
                                                  bonf, alpha,
                                                  max_points=report_max_points,
                                                  rasterize_scatter=report_rasterize_scatter,
                                                  random_seed=report_seed)
                    if fig:
                        pdf.savefig(fig, dpi=report_pdf_dpi)
                        plt.close(fig)

            # --- QQ per model ---
            for method in model_methods:
                p_col = f'p_value_{method}'
                if p_col in results_df.columns:
                    pvals = results_df[p_col].dropna().values
                    fig = _plot_qq_for_pdf(
                        pvals,
                        f'QQ Plot – {method.upper()}',
                        max_points=report_max_points,
                        rasterize_scatter=report_rasterize_scatter,
                        random_seed=report_seed,
                    )
                    if fig:
                        pdf.savefig(fig, dpi=report_pdf_dpi)
                        plt.close(fig)

            # --- Comparative QQ ---
            fig = _plot_comparative_qq_for_pdf(
                qq_data,
                max_points=report_max_points,
                rasterize_scatter=report_rasterize_scatter,
                random_seed=report_seed,
            )
            if fig:
                pdf.savefig(fig, dpi=report_pdf_dpi)
                plt.close(fig)

            # --- Volcano plot ---
            p_col_v = f'p_value_{target_model}'
            b_col_v = f'beta_{target_model}'
            if p_col_v in results_df.columns and b_col_v in results_df.columns:
                fig = _plot_volcano_for_pdf(
                    results_df,
                    p_col_v,
                    b_col_v,
                    f'Volcano Plot – {target_model.upper()}',
                    max_points=report_max_points,
                    rasterize_scatter=report_rasterize_scatter,
                    random_seed=report_seed,
                )
                if fig:
                    pdf.savefig(fig, dpi=report_pdf_dpi)
                    plt.close(fig)

            # --- Top Loci Table ---
            _add_top_loci_page(pdf, top_loci, target_model)

            # --- Pathway Enrichment ---
            _add_enrichment_page(pdf, enrichment_df)

            # --- COJO Conditional Analysis ---
            _add_cojo_page(pdf, cojo_results)

            # --- Polygenic Risk Score ---
            _add_prs_page(pdf, prs_summary)

            # --- GWAS Catalog Cross-Reference ---
            _add_gwas_catalog_page(pdf, catalog_results)

            # --- Session Info ---
            _add_session_page(pdf, config,
                              elapsed_seconds=report_data.get('elapsed'),
                              peak_memory=report_data.get('peak_memory'))

    finally:
        # Restore original backend so notebook inline plots keep working
        try:
            matplotlib.rcParams['pdf.compression'] = _orig_pdf_comp
            plt.switch_backend(orig_backend)
        except Exception:
            pass

    # Also write a brief markdown summary for backward compatibility
    if output_path:
        try:
            n_samples = config.get('n_samples', '?')
            n_snps = config.get('n_snps', '?')
            with open(output_path, 'w') as f:
                f.write(f"# GWAS Summary Report\n\n")
                f.write(f"**Trait:** {config.get('TRAIT_NAME', '?')}  \n")
                f.write(f"**Species:** {config.get('SPECIES', '?')}  \n")
                f.write(f"**Trait Type:** {config.get('TRAIT_TYPE', '?')}  \n")
                f.write(f"**Samples:** {n_samples}  |  **SNPs (post-QC):** {n_snps}  \n")
                f.write(f"**Date:** {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}  \n\n")
                if best_model:
                    lam = lambda_values.get(best_model)
                    lam_s = f' (λ={lam:.3f})' if lam else ''
                    f.write(f"**Recommended model:** {best_model}{lam_s}  \n\n")
                if summary_df is not None and not summary_df.empty:
                    f.write("## Model Comparison\n\n")
                    f.write(summary_df.to_markdown(index=False))
                    f.write("\n\n")
                if top_loci is not None and not top_loci.empty:
                    f.write(f"## Top Loci ({target_model.upper()}, {len(top_loci)} loci)\n\n")
                    show = top_loci.head(15)
                    cols = [c for c in ['SNP_ID', 'chrom', 'pos', f'p_value_{target_model}',
                                        'Gene_Symbol'] if c in show.columns]
                    if cols:
                        f.write(show[cols].to_markdown(index=False))
                    f.write("\n\n")
                f.write(f"\n---\n*Full PDF report:* `{os.path.basename(pdf_path)}`\n")
        except Exception as e:
            print(f"Markdown summary warning: {e}")

    sz = os.path.getsize(pdf_path) if os.path.exists(pdf_path) else 0
    sz_str = f'{sz / 1024:.1f} KB' if sz < 1024 ** 2 else f'{sz / 1024 ** 2:.1f} MB'
    print(f"PDF report saved: {pdf_path} ({sz_str})")
    return pdf_path
