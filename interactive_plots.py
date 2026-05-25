import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from scipy import stats
import logging
import seaborn as sns
import matplotlib.pyplot as plt

# Try to import nbformat for Plotly static image export validation
try:
    import nbformat
except ImportError:
    logging.warning("nbformat not installed. Interactive plots might not render correctly in some environments.")


def create_interactive_manhattan(df, p_col, chrom_col='chrom', pos_col='pos', snp_col='SNP_ID', 
                               title="Manhattan Plot", suggestive_threshold=1e-5, bonferroni_threshold=None):
    """
    Creates an interactive Manhattan plot using Plotly.
    """
    df_plot = df.dropna(subset=[p_col, chrom_col, pos_col]).copy()
    df_plot['-log10p'] = -np.log10(df_plot[p_col].astype(float) + 1e-300)
    
    # Ensure chromosome is string and sorted
    df_plot[chrom_col] = df_plot[chrom_col].astype(str)
    
    # Create a continuous X axis for the plot
    # We need to calculate cumulative offsets for each chromosome
    chrom_sizes = df_plot.groupby(chrom_col)[pos_col].max()
    
    # Sort chromosomes naturally
    try:
        from natsort import natsorted
        sorted_chroms = natsorted(df_plot[chrom_col].unique())
    except ImportError:
        sorted_chroms = sorted(df_plot[chrom_col].unique())
        
    chrom_offsets = {}
    current_offset = 0
    tick_positions = []
    tick_labels = []
    
    for chrom in sorted_chroms:
        chrom_offsets[chrom] = current_offset
        size = chrom_sizes[chrom]
        tick_positions.append(current_offset + size/2)
        tick_labels.append(chrom)
        current_offset += size
        
    df_plot['global_pos'] = df_plot.apply(lambda row: chrom_offsets[row[chrom_col]] + row[pos_col], axis=1)
    
    # Downsample non-significant SNPs for performance if dataset is huge
    # Keep all significant SNPs, sample 10% of others
    sig_mask = df_plot['-log10p'] > -np.log10(suggestive_threshold)
    non_sig_df = df_plot[~sig_mask]
    if len(non_sig_df) > 50000:
        non_sig_df = non_sig_df.sample(n=50000, random_state=42)
    
    df_reduced = pd.concat([df_plot[sig_mask], non_sig_df])
    
    fig = px.scatter(
        df_reduced, 
        x='global_pos', 
        y='-log10p',
        color=chrom_col,
        hover_data=[snp_col, chrom_col, pos_col, p_col],
        title=title,
        labels={'global_pos': 'Chromosome', '-log10p': '-log10(p-value)'},
        category_orders={chrom_col: sorted_chroms}
    )
    
    # Update x-axis to show chromosome names
    fig.update_layout(
        xaxis=dict(
            tickmode='array',
            tickvals=tick_positions,
            ticktext=tick_labels,
            showgrid=False
        ),
        showlegend=True
    )
    
    # Add threshold lines
    if suggestive_threshold:
        fig.add_hline(y=-np.log10(suggestive_threshold), line_dash="dot", 
                     annotation_text=f"Suggestive ({suggestive_threshold:.1e})", 
                     line_color="blue")
        
    if bonferroni_threshold:
        fig.add_hline(y=-np.log10(bonferroni_threshold), line_dash="dash", 
                     annotation_text=f"Bonferroni ({bonferroni_threshold:.1e})", 
                     line_color="red")
        
    return fig

def create_interactive_qq(p_values, title="QQ Plot"):
    """
    Creates an interactive QQ plot with confidence intervals using Plotly.
    """
    p_clean = pd.Series(p_values).dropna()
    if p_clean.empty:
        return None
        
    n = len(p_clean)
    observed = -np.log10(np.sort(p_clean))
    expected = -np.log10(np.arange(1, n + 1) / (n + 1))
    
    # Calculate Lambda GC
    chisq = stats.chi2.ppf(1 - p_clean, 1)
    lambda_gc = np.median(chisq) / stats.chi2.ppf(0.5, 1)
    
    # Confidence Intervals (95%)
    # Beta distribution quantiles
    k = np.arange(1, n + 1)
    beta_lower = stats.beta.ppf(0.025, k, n - k + 1)
    beta_upper = stats.beta.ppf(0.975, k, n - k + 1)
    
    ci_lower = -np.log10(beta_lower)
    ci_upper = -np.log10(beta_upper)
    
    # Downsample for plotting if too large
    if n > 10000:
        indices = np.unique(np.floor(np.logspace(0, np.log10(n-1), 10000)).astype(int))
        observed = observed[indices]
        expected = expected[indices]
        ci_lower = ci_lower[indices]
        ci_upper = ci_upper[indices]
    
    fig = go.Figure()
    
    # Confidence Interval (Shaded Area)
    fig.add_trace(go.Scatter(
        x=np.concatenate([expected, expected[::-1]]),
        y=np.concatenate([ci_upper, ci_lower[::-1]]),
        fill='toself',
        fillcolor='rgba(200, 200, 200, 0.5)',
        line=dict(color='rgba(255,255,255,0)'),
        hoverinfo="skip",
        name='95% CI'
    ))
    
    # Observed vs Expected
    fig.add_trace(go.Scatter(
        x=expected, 
        y=observed, 
        mode='markers',
        marker=dict(size=5, color='blue'),
        name='Observed'
    ))
    
    # Diagonal line
    max_val = max(expected.max(), observed.max())
    fig.add_trace(go.Scatter(
        x=[0, max_val],
        y=[0, max_val],
        mode='lines',
        line=dict(color='red', dash='dash'),
        name='Expected'
    ))
    
    fig.update_layout(
        title=f"{title} (λ_GC = {lambda_gc:.3f})",
        xaxis_title="Expected -log10(p)",
        yaxis_title="Observed -log10(p)",
        hovermode="closest"
    )
    
    return fig

def perform_distance_clumping(df, p_col, chrom_col='chrom', pos_col='pos', snp_col='SNP_ID', 
                            window_kb=250, p_threshold=1e-5):
    """
    Performs distance-based clumping to identify independent significant loci.
    """
    sig_df = df[df[p_col] < p_threshold].copy()
    
    if sig_df.empty:
        return pd.DataFrame()
        
    # Sort by p-value (most significant first)
    sig_df = sig_df.sort_values(p_col)
    
    clumped_snps = []
    
    while not sig_df.empty:
        # Take the most significant SNP (index 0 after sort)
        lead_snp = sig_df.iloc[0]
        clumped_snps.append(lead_snp)
        
        # Remove SNPs within window on the same chromosome
        mask = (
            (sig_df[chrom_col] == lead_snp[chrom_col]) & 
            (abs(sig_df[pos_col] - lead_snp[pos_col]) <= window_kb * 1000)
        )
        
        # Drop these SNPs from the pool
        sig_df = sig_df[~mask]
        
    return pd.DataFrame(clumped_snps)

def create_interactive_volcano(df, p_col, effect_col='beta', snp_col='SNP_ID', 
                             title="Volcano Plot", p_threshold=1e-5):
    """
    Creates an interactive Volcano plot (Effect Size vs -log10(p-value)).
    """
    df_plot = df.dropna(subset=[p_col, effect_col]).copy()
    df_plot['-log10p'] = -np.log10(df_plot[p_col].astype(float) + 1e-300)
    
    # Color coding
    df_plot['color'] = 'Not Significant'
    df_plot.loc[df_plot['-log10p'] > -np.log10(p_threshold), 'color'] = 'Significant'
    
    fig = px.scatter(
        df_plot, 
        x=effect_col, 
        y='-log10p',
        color='color',
        hover_data=[snp_col, p_col, effect_col],
        title=title,
        labels={effect_col: 'Effect Size (Beta)', '-log10p': '-log10(p-value)'},
        color_discrete_map={'Significant': 'red', 'Not Significant': 'grey'}
    )
    
    fig.add_hline(y=-np.log10(p_threshold), line_dash="dash", line_color="blue", annotation_text=f"Threshold ({p_threshold})")
    
    return fig

def create_ld_heatmap(geno_df, snp_list, title="LD Heatmap"):
    """
    Creates a Linkage Disequilibrium (LD) heatmap for a list of SNPs.
    Uses correlation (r^2) between genotypes.
    """
    if len(snp_list) < 2:
        return None
        
    # Subset genotype data
    # Ensure SNPs exist in the dataframe
    valid_snps = [snp for snp in snp_list if snp in geno_df.columns]
    
    if len(valid_snps) < 2:
        return None
        
    sub_geno = geno_df[valid_snps].astype(float)
    
    # Calculate correlation matrix (r)
    corr_matrix = sub_geno.corr()
    
    # Calculate r^2
    ld_matrix = corr_matrix ** 2
    
    # Create heatmap using Plotly for interactivity
    fig = go.Figure(data=go.Heatmap(
        z=ld_matrix.values,
        x=ld_matrix.columns,
        y=ld_matrix.index,
        colorscale='Viridis',
        zmin=0, zmax=1
    ))
    
    fig.update_layout(
        title=title,
        xaxis_title="SNP",
        yaxis_title="SNP",
        width=800,
        height=800
    )
    
    return fig


# ======================================================================
# LD-BASED CLUMPING (r²)
# ======================================================================

def perform_ld_clumping(results_df, G_qc, snp_ids, p_col, r2_threshold=0.1,
                        window_kb=500, p_threshold=1e-5,
                        chrom_col='chrom', pos_col='pos', snp_col='SNP_ID'):
    """
    LD-based clumping: keeps independent lead SNPs where no other significant
    SNP within *window_kb* has r² > r2_threshold with the lead.

    Parameters
    ----------
    results_df  : DataFrame with p-values, chrom, pos
    G_qc        : numpy array (n_samples x n_snps) -- QC'd genotypes
    snp_ids     : array of SNP IDs matching G_qc columns
    p_col       : column name for p-values
    r2_threshold: max r² allowed between lead and clump members (default 0.1)
    window_kb   : clumping window in kilobases (default 500)
    p_threshold : significance threshold for candidate SNPs

    Returns
    -------
    DataFrame of clumped lead SNPs.
    """
    df = results_df.dropna(subset=[p_col]).copy()
    sig = df[df[p_col] < p_threshold].sort_values(p_col)
    if sig.empty:
        return pd.DataFrame()

    # Build SNP-ID to column-index map
    id_to_idx = {sid: i for i, sid in enumerate(snp_ids)}

    clumped = []
    removed = set()
    window_bp = window_kb * 1000

    for _, lead in sig.iterrows():
        lead_id = lead[snp_col]
        if lead_id in removed:
            continue
        if lead_id not in id_to_idx:
            continue
        clumped.append(lead)
        lead_idx = id_to_idx[lead_id]
        lead_chr = str(lead[chrom_col])
        lead_pos = int(lead[pos_col])
        lead_g   = G_qc[:, lead_idx].astype(float)
        lead_g_filled = np.where(np.isnan(lead_g), 0, lead_g)
        lg_std = lead_g_filled.std()
        if lg_std == 0:
            continue

        # Find candidates on same chromosome within window
        cand_mask = (
            (sig[chrom_col].astype(str) == lead_chr) &
            (abs(sig[pos_col].astype(int) - lead_pos) <= window_bp) &
            (~sig[snp_col].isin(removed)) &
            (sig[snp_col] != lead_id)
        )
        for _, cand in sig[cand_mask].iterrows():
            cand_id = cand[snp_col]
            if cand_id in removed or cand_id not in id_to_idx:
                continue
            cand_idx = id_to_idx[cand_id]
            cand_g = G_qc[:, cand_idx].astype(float)
            cand_g_filled = np.where(np.isnan(cand_g), 0, cand_g)
            cg_std = cand_g_filled.std()
            if cg_std == 0:
                removed.add(cand_id)
                continue
            r = np.corrcoef(lead_g_filled, cand_g_filled)[0, 1]
            if r ** 2 > r2_threshold:
                removed.add(cand_id)

    return pd.DataFrame(clumped).reset_index(drop=True) if clumped else pd.DataFrame()


# ======================================================================
# COJO -- Conditional & Joint analysis (simplified, pure-Python)
# ======================================================================

def cojo_conditional(results_df, G_qc, y, snp_ids, p_col,
                     lead_snp_ids, window_kb=1000, p_cojo=1e-5,
                     max_secondary=5, covariates=None):
    """
    Simple conditional analysis -- for each locus, condition on the lead SNP
    and re-test nearby SNPs to find secondary signals.

    Parameters
    ----------
    results_df   : DataFrame with p-values
    G_qc         : genotype array
    y            : phenotype array
    snp_ids      : SNP ID array
    p_col        : p-value column
    lead_snp_ids : list of lead SNP IDs from clumping
    window_kb    : window around each lead to re-test
    p_cojo       : threshold for declaring secondary signals
    max_secondary: max secondary hits per locus
    covariates   : optional covariate array

    Returns
    -------
    DataFrame with columns [Locus_Lead, Secondary_SNP, chrom, pos, p_conditional, beta_cond, se_cond]
    """
    from scipy import stats as _stats
    id_to_idx = {sid: i for i, sid in enumerate(snp_ids)}
    records = []

    for lead_id in lead_snp_ids:
        if lead_id not in id_to_idx:
            continue
        lead_idx = id_to_idx[lead_id]
        lead_row = results_df[results_df['SNP_ID'] == lead_id]
        if lead_row.empty:
            continue
        lead_chr = str(lead_row['chrom'].iloc[0])
        lead_pos = int(lead_row['pos'].iloc[0])
        window_bp = window_kb * 1000

        # Nearby SNPs on same chromosome
        nearby = results_df[
            (results_df['chrom'].astype(str) == lead_chr) &
            (abs(results_df['pos'].astype(int) - lead_pos) <= window_bp) &
            (results_df['SNP_ID'] != lead_id)
        ]
        if nearby.empty:
            continue

        # Build base design: intercept + lead SNP genotype [+ covariates]
        g_lead = G_qc[:, lead_idx].astype(float)
        parts = [np.ones((len(y), 1)), g_lead.reshape(-1, 1)]
        if covariates is not None:
            parts.append(covariates)
        X_base = np.column_stack(parts)

        n_found = 0
        for _, row in nearby.iterrows():
            if n_found >= max_secondary:
                break
            sid = row['SNP_ID']
            if sid not in id_to_idx:
                continue
            sidx = id_to_idx[sid]
            g_snp = G_qc[:, sidx].astype(float)
            X = np.column_stack([X_base, g_snp])
            try:
                import statsmodels.api as sm
                m = sm.OLS(y, X).fit()
                pv = m.pvalues[-1]
                if pv < p_cojo:
                    records.append({
                        'Locus_Lead': lead_id,
                        'Secondary_SNP': sid,
                        'chrom': row['chrom'],
                        'pos': row['pos'],
                        'p_conditional': pv,
                        'beta_cond': m.params[-1],
                        'se_cond': m.bse[-1],
                    })
                    n_found += 1
            except:
                pass

    return pd.DataFrame(records) if records else pd.DataFrame()


# ======================================================================
# LOCUS-ZOOM PLOT
# ======================================================================

def create_locus_zoom(results_df, G_qc, snp_ids, lead_snp_id, p_col,
                      window_kb=500, chrom_col='chrom', pos_col='pos',
                      snp_col='SNP_ID', gene_annotations=None):
    """
    Interactive regional association plot (locus-zoom style).
    Color SNPs by LD (r²) with the lead variant.

    Parameters
    ----------
    results_df     : DataFrame with p-values, chrom, pos
    G_qc           : genotype matrix
    snp_ids        : SNP ID array
    lead_snp_id    : ID of the lead SNP
    p_col          : p-value column
    window_kb      : plot window in kb
    gene_annotations : optional DataFrame with gene_start, gene_end, gene_name for overlay

    Returns
    -------
    plotly Figure
    """
    id_to_idx = {sid: i for i, sid in enumerate(snp_ids)}
    lead_row = results_df[results_df[snp_col] == lead_snp_id]
    if lead_row.empty or lead_snp_id not in id_to_idx:
        return None

    lead_chr = str(lead_row[chrom_col].iloc[0])
    lead_pos = int(lead_row[pos_col].iloc[0])
    lead_idx = id_to_idx[lead_snp_id]
    window_bp = window_kb * 1000

    # Extract regional SNPs
    region = results_df[
        (results_df[chrom_col].astype(str) == lead_chr) &
        (abs(results_df[pos_col].astype(int) - lead_pos) <= window_bp) &
        results_df[p_col].notna()
    ].copy()

    if region.empty:
        return None

    # Compute r² with lead SNP
    lead_g = G_qc[:, lead_idx].astype(float)
    lead_g_filled = np.where(np.isnan(lead_g), 0, lead_g)
    lg_std = lead_g_filled.std()

    r2_vals = []
    for _, row in region.iterrows():
        sid = row[snp_col]
        if sid == lead_snp_id:
            r2_vals.append(1.0)
        elif sid in id_to_idx and lg_std > 0:
            cg = G_qc[:, id_to_idx[sid]].astype(float)
            cg = np.where(np.isnan(cg), 0, cg)
            if cg.std() > 0:
                r = np.corrcoef(lead_g_filled, cg)[0, 1]
                r2_vals.append(r ** 2)
            else:
                r2_vals.append(0.0)
        else:
            r2_vals.append(0.0)

    region['r2_with_lead'] = r2_vals
    region['-log10p'] = -np.log10(region[p_col].astype(float) + 1e-300)

    # LD colour bins
    def _ld_color(r2):
        if r2 >= 0.8: return 'red'
        if r2 >= 0.6: return 'orange'
        if r2 >= 0.4: return 'green'
        if r2 >= 0.2: return 'skyblue'
        return 'navy'

    region['ld_color'] = region['r2_with_lead'].apply(_ld_color)
    region['ld_bin'] = pd.cut(region['r2_with_lead'],
                              bins=[-0.01, 0.2, 0.4, 0.6, 0.8, 1.01],
                              labels=['0-0.2', '0.2-0.4', '0.4-0.6', '0.6-0.8', '0.8-1.0'])

    color_map = {'0-0.2': 'navy', '0.2-0.4': 'skyblue', '0.4-0.6': 'green',
                 '0.6-0.8': 'orange', '0.8-1.0': 'red'}

    fig = px.scatter(
        region, x=pos_col, y='-log10p', color='ld_bin',
        color_discrete_map=color_map,
        hover_data=[snp_col, p_col, 'r2_with_lead'],
        title=f"Locus Zoom: {lead_snp_id} (chr{lead_chr})",
        labels={pos_col: f'Position on chr{lead_chr} (bp)', '-log10p': '-log10(p)'},
        category_orders={'ld_bin': ['0-0.2', '0.2-0.4', '0.4-0.6', '0.6-0.8', '0.8-1.0']},
    )

    # Mark lead SNP with a diamond
    lead_data = region[region[snp_col] == lead_snp_id]
    if not lead_data.empty:
        fig.add_trace(go.Scatter(
            x=lead_data[pos_col], y=lead_data['-log10p'],
            mode='markers', marker=dict(size=14, symbol='diamond', color='purple', line=dict(width=2, color='black')),
            name='Lead SNP', showlegend=True,
            hovertext=lead_snp_id,
        ))

    # Gene annotations if available
    if gene_annotations is not None and not gene_annotations.empty:
        y_gene = region['-log10p'].max() * -0.05
        for _, gene in gene_annotations.iterrows():
            fig.add_shape(type='line', x0=gene.get('gene_start', 0), x1=gene.get('gene_end', 0),
                          y0=y_gene, y1=y_gene, line=dict(color='darkgreen', width=4))
            fig.add_annotation(x=(gene.get('gene_start',0)+gene.get('gene_end',0))/2,
                             y=y_gene - 0.3, text=gene.get('gene_name',''), showarrow=False,
                             font=dict(size=8))

    fig.update_layout(legend_title='r² with lead', template='plotly_white')
    return fig


# ======================================================================
# GWAS CATALOG LOOKUP (REST API)
# ======================================================================

def gwas_catalog_lookup(top_loci_df, window_bp=50000, chrom_col='chrom', pos_col='pos',
                        snp_col='SNP_ID', max_hits=5):
    """
    Query the NHGRI-EBI GWAS Catalog REST API for known associations near top hits.

    Parameters
    ----------
    top_loci_df : DataFrame with chrom and pos columns
    window_bp   : search window in basepairs around each locus
    max_hits    : max catalog results to keep per locus

    Returns
    -------
    DataFrame with columns [Query_SNP, chrom, pos, Catalog_RSid, Catalog_Trait, Catalog_P, Distance_bp]
    """
    import requests as _req, time as _t

    records = []
    base = "https://www.ebi.ac.uk/gwas/rest/api"

    for _, row in top_loci_df.iterrows():
        chrom = str(row[chrom_col]).replace('chr', '')
        pos = int(row[pos_col])
        start = max(0, pos - window_bp)
        end = pos + window_bp
        query_snp = row[snp_col]

        url = f"{base}/associations/search/findByChromosomeAndBpRange?chromosome={chrom}&start={start}&end={end}"
        try:
            r = _req.get(url, headers={"Accept": "application/json"}, timeout=10)
            if r.ok:
                data = r.json()
                assocs = data.get('_embedded', {}).get('associations', [])
                for a in assocs[:max_hits]:
                    trait_name = '; '.join([t.get('trait', '') for t in a.get('efoTraits', [])])
                    pv = a.get('pvalue', None)
                    strongest = ''
                    for sra in a.get('strongestRiskAlleles', []):
                        strongest = sra.get('riskAlleleName', '')
                        break
                    mid = (start + end) // 2
                    records.append({
                        'Query_SNP': query_snp,
                        'chrom': chrom, 'pos': pos,
                        'Catalog_RSid': strongest,
                        'Catalog_Trait': trait_name,
                        'Catalog_P': pv,
                        'Distance_bp': abs(pos - mid),
                    })
            _t.sleep(0.2)  # rate limit
        except Exception:
            pass

    return pd.DataFrame(records) if records else pd.DataFrame()
