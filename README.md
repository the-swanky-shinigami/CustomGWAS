# CustomGWAS: An Interactive End-to-End GWAS Pipeline in Python

[![Python Version](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

CustomGWAS is a comprehensive, interactive, and user-friendly Jupyter Notebook for conducting a complete Genome-Wide Association Study (GWAS). It transforms the traditionally complex, multi-tool "black-box" process of GWAS into a single, transparent "white-box" environment.

This pipeline is designed to democratize GWAS, enabling researchers of all computational skill levels to perform a robust, publication-quality analysis—from raw data to final plots—with confidence and clarity.

##  Key Features

* ✅ **End-to-End Workflow**: A single notebook manages all steps: data loading, quality control, association testing, and visualization.
* 🧠 **Automated Data Handling**: Intelligently handles both CSV and PLINK formats. Automatically converts large CSV files to the memory-efficient PLINK binary format to prevent crashes and optimize performance.
* 📊 **Interactive Quality Control (QC)**: Implements a rigorous QC pipeline with interactive widgets that provide detailed reports on sample and SNP dropout at each step.
* 📈 **Enhanced Phenotype Analysis**: Includes checks for phenotype distribution, automated data transformation for skewed data, and outlier removal to ensure robust results.
* 🚀 **Multi-Model Framework**: Runs up to four different association models in parallel to provide a comprehensive view of the results:
    * **Naive GLM**: A basic model for baseline comparison.
    * **PCA-Corrected GLM**: Corrects for population stratification.
    * **Linear Mixed Model (LMM)**: The industry standard, correcting for both population structure and cryptic relatedness.
    * **Leave-One-Chromosome-Out (LOCO) LMM**: The gold standard, avoiding proximal contamination for more accurate results.
* 🎨 **Publication-Ready Visuals**: Automatically generates high-quality Manhattan plots and QQ plots for each model.
* reproducibility_and_transparency.png: 

## Installation

To run the CustomGWAS pipeline, you will need Python 3.8+ and the required packages. A C++ compiler is also required for some dependencies.

### 1. Clone the Repository

```bash
git clone [https://github.com/the-swanky-shinigami/CustomGWAS.git](https://github.com/the-swanky-shinigami/CustomGWAS.git)
cd CustomGWAS
