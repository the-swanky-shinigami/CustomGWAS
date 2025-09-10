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
* ⚙️ **Fully Customizable** : Customize the notebook according to your needs and manipulate the parameters in the beginning to achieve desired results
 

## Installation

To run the CustomGWAS pipeline, you will need Python 3.8+ and the required packages. A C++ compiler is also required for some dependencies.

### 1. Clone the Repository

```bash
git clone [https://github.com/the-swanky-shinigami/CustomGWAS.git](https://github.com/the-swanky-shinigami/CustomGWAS.git)
cd CustomGWAS
```

### 2. Install Dependencies

It is highly recommended to use a conda environment to manage the dependencies.

```bash
# Create and activate a new conda environment
conda create -n gwas_env python=3.8
conda activate gwas_env

# Install dependencies from requirements.txt
pip install -r requirements.txt
```
### 3. Install PLINK (Required for large CSV files)

CustomGWAS leverages PLINK for efficient, out-of-core handling of large genotype files in CSV format. You must install PLINK and ensure it is available in your system's PATH.
* **Linux/macOS**:
   ```bash
   # Using conda
   conda install -c bioconda plink
   ```
* **Windows**: Download the appropriate PLINK 1.9 binary from the official website and add its location to your system's PATH environment variable.

## How to Use the Pipeline

1. **Launch Jupyter**: Make sure your gwas_env conda environment is active and start Jupyter Lab or Jupyter Notebook:
   ```bash
   jupyter lab
   ```
2. **Open the Notebook**: Open the ```CUSTOM_GWAS_NOTEBOOK.ipynb``` file.
3. **Configure Your Analysis**: Go to the "**1. SETUP AND CONFIGURATION**" cell. This is the only cell you need to edit.
   * Set the file paths for your genotype and phenotype data. You can use the provided ```Sample_Genotype_Data.csv``` and ```Sample_Phenotype_Data.csv``` to get started immediately.
   * Adjust QC thresholds (e.g., MAF, HWE, missingness) as needed.
   * Select which GWAS models you want to run (```True```/```False```).
   * Set your output directory and project name.
4. **Run the Pipeline**: Click "```Run```" > "```Run All Cells```" in the Jupyter menu. The notebook will execute the entire workflow from start to finish.
5. **Review Your Results**: The notebook will display interactive QC tables, generate Manhattan and QQ plots, and save all results (including top SNP lists) to your specified output directory.

## Included Sample Data

To help you get started, this repository includes sample datasets:

   * ```Sample_Genotype_Data.csv```: A sample genotype file containing SNP data for a number of individuals.

   * ```Sample_Phenotype_Data.csv```: A corresponding phenotype file with trait values for the same individuals.

These files can be configured in the notebook, so you can run the pipeline out-of-the-box to see how it works.

## How to Cite

If you use CustomGWAS in your research, please cite our manuscript:

[Details will be available once our work is published, for the time being you can use the repo link to cite our work]

## License
This project is licensed under the MIT License - see the ```LICENSE.md``` file for details.
