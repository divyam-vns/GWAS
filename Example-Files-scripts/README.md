# GWAS Pipeline

This is a GitHub-ready pipeline for performing Genome-Wide Association Studies (GWAS) from raw genotype and phenotype data.

## Requirements
- PLINK (v1.9 or v2.0)
- Python 3.10+
- Python packages: see `requirements.txt`

## Pipeline Steps
* 1. **Quality Control**
   ```bash
     bash scripts/qc_plink.sh
* 2. **PCA for population stratification**
  ```bash
     bash scripts/pca_plink.sh
* 3. **GWAS Analysis & Visualization**
    ```bash
    python scripts/gwas_analysis.py
## Output
```results/gwas_results.csv```: GWAS results with p-values and FDR

```results/manhattan_plot.png```: Manhattan plot

```results/qq_plot.png```: QQ plot
