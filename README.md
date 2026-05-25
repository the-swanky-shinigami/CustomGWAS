# CustomGWAS Pipeline v3

CustomGWAS Pipeline v3 is an end-to-end GWAS (Genome-Wide Association Study) workflow implemented as an interactive Python Jupyter notebook. It provides scalable genotype ingestion, robust quality control, multi-model association testing, post-GWAS follow-up, and reproducible reporting.

This pipeline wraps high-performance libraries (like `pandas-plink`) and mathematically rigorous linear models in a user-friendly, transparent environment, effectively abstracting away the typical boilerplate code of GWAS while allowing researchers complete control over their parameters.

## Core Files
- `GWAS_Pipeline_v3.ipynb`: The primary executable notebook.
- `gwas_analytics.py`: Generates the automated PDF and Markdown summary reports.
- `interactive_plots.py`: Renders interactive Plotly graphs (Volcano plots, Manhattan plots, LD Heatmaps) saved directly as HTML artifacts.
- `r_blink_wrapper.py` & `run_blink.R`: External bridge scripts that safely translate massive python-loaded arrays to R environments for execution using `BLINK`, `FarmCPU`, and `mrMLM`.

---

## 🚀 Feature Coverage

### 1. Genotype Ingestion & Preprocessing
- **Multi-format Support**: Seamlessly ingest CSV, `.bed/.bim/.fam` PLINK triplets, or VCF genotypes.
- **Auto-conversion**: Includes an ultra-optimized CSV-to-PLINK conversion block for massive raw text datasets.
- **Alignment**: Automatically intersects and aligns genotype samples, phenotype records, and user-provided covariates.

### 2. Quality Control (QC)
- Sample and SNP call-rate filtering.
- Heterozygosity outlier filtering.
- Minor Allele Frequency (MAF) filtering.
- Hardy-Weinberg Equilibrium (HWE) testing and filtering.
- **Imputation**: Support for Mean, Median, KNN, or external Beagle imputation.
- **Phenotype transformation**: Detects skewness and gracefully applies log/Box-Cox transformations if enabled.

### 3. Population Structure & Kinship Correction
- **PCA**: Automatic LD-pruning and calculation of Principal Components to control for population stratification.
- **Kinship Matrix**: VanRaden-based genomic relationship matrix calculation for mixed models.

### 4. Association Testing Models
The pipeline calculates associations across an array of nested linear models and multi-locus approaches:
- **Naive GLM** (baseline)
- **PCA-corrected GLM** (controls for population structure)
- **EMMAX / LMM** (controls for both population structure and cryptic relatedness/kinship)
- **LOCO-LMM** (Leave-One-Chromosome-Out LMM to prevent proximal contamination)
- **Stepwise MLM** 
- **Multi-locus Models**: Native integration with R packages `BLINK`, `FarmCPU`, and `mrMLM`.

### 5. Post-GWAS Analytics
- **Model Recommendation**: Evaluates the Genomic Inflation Factor ($\lambda$) for all successfully completed models and recommends the model whose $\lambda$ is closest to $1.0$ (best calibration).
- **False Discovery Control**: Bonferroni and Benjamini-Hochberg FDR ($q$-value) thresholds.
- **Clumping**: LD-based clumping to collapse thousands of significant SNPs into independent Top Loci.
- **Annotation**: Automatic queries against biological databases (e.g. MyGene.info) to annotate nearest genes for significant loci.
- **Polygenic Risk Scoring (PRS)**: Clump-and-threshold PRS evaluation with automatic holdout scoring.

---

## 🛠️ Installation & Usage

### Option 1: Docker (Highly Recommended)
Docker abstracts away the underlying operating system. Running this pipeline in Docker guarantees 100% feature parity across **Windows**, **macOS**, and **Linux** systems without manually managing R versions, Python environments, Java installations, or PLINK binaries.

1. Ensure [Docker Desktop](https://www.docker.com/products/docker-desktop) is installed and running on your system.
2. Clone or transfer the directory to your target machine.
3. Open a terminal in the directory and run:
   ```bash
   # Run this ONLY the first time or if you change the Dockerfile/requirements.txt
   docker compose build 
   
   # Run this to start the environment (you can just run this next time you open it)
   docker compose up
   ```
4. Open the Jupyter link provided in the terminal (e.g., `http://127.0.0.1:8888/...`).

**Caveat - Docker Memory Constraints:** 
For exceptionally large datasets (e.g. Rice 3KRG, containing over 3000 samples and millions of SNPs), Docker's default memory limits (usually 2-8 GB) will trigger `OOM (Out Of Memory)` Kernel crashes or `DLASCL` math errors. 
- **Solution**: Open Docker Desktop Settings -> Resources, and increase your Memory Limit to at least 16GB (or 32GB+ for massive arrays), and increase Swap. 

### Option 2: Native Host (Python Virtual Environment)
If you prefer running natively, you must have Python 3.10+, R, Java, and PLINK installed on your host OS.

1. Install system dependencies:
   - [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/) (must be accessible in `$PATH` or in the working directory).
   - [Java](https://www.java.com/en/) (required if Beagle imputation is enabled).
   - [R](https://www.r-project.org/) and `Rscript` (required for BLINK, FarmCPU, mrMLM).
2. Install Python dependencies:
   ```bash
   python3 -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate
   pip install -r requirements.txt
   ```
3. Install R dependencies:
   ```R
   install.packages(c("MASS", "mrMLM", "BLINK", "FarmCPUpp"))
   ```
4. Run `jupyter notebook` and open `GWAS_Pipeline_v3.ipynb`.

---

## 📁 Transferring the Pipeline
If you wish to move this tool to another computer and run it via Docker, simply copy the entire root directory, or at a minimum, ensure the following files are preserved:

**Required Files for Docker Execution**:
- `Dockerfile`
- `docker-compose.yml`
- `requirements.txt`
- `docker/install_r_packages.R`
- `docker/entrypoint.sh`

**Required Files for Pipeline Execution**:
- `GWAS_Pipeline_v3.ipynb`
- `gwas_analytics.py`
- `interactive_plots.py`
- `r_blink_wrapper.py` & `run_blink.R`
- `beagle.27Feb25.75f.jar` (if imputing)
- Your `runs/` directory (if you wish to keep past logs/results)
- Your input datasets (`.csv`, `.bed/.bim/.fam`)

*(Note: There are absolutely no OS-specific complications when using Docker. The pipeline will operate identically on Windows, Linux, and macOS. The wrapper scripts intelligently map execution paths automatically).*

---

## 📝 Outputs and Reporting

All outputs for a single execution are isolated within timestamped directories under `runs/` (e.g. `runs/Anandan_GWAS_v3_output_20260428_112911`). 

Inside, you will find:
1. `*_gwas_results.csv`: The master table containing all calculated p-values and effects for all models.
2. `*_top_loci.csv`: The LD-clumped significant hits across all chromosomes.
3. `*_Summary_Report.md` / `*.pdf`: Auto-generated, presentation-ready documents outlining the quality control steps, genomic inflation, and final findings.
4. `*.html`: Interactive plots.
5. `pipeline.log`: Verbose trace of the entire execution.
