# Overview of GWAS Pipeline

## 1. A typical GWAS pipeline involves these steps:

### 1. Data preparation

* Genotype data (SNPs) → PLINK format (.bed/.bim/.fam)
* Phenotype data (traits) → CSV/TSV
* Covariates (age, sex, PCs, batch effects) → CSV/TSV

### 2. Quality Control (QC)

* Remove SNPs with low call rate
* Remove individuals with high missingness
* Remove SNPs failing Hardy-Weinberg equilibrium
* Filter by minor allele frequency (MAF)

### 3. Population Stratification

* PCA on genotype data to correct for population structure

### 4. Association Testing

* Linear or logistic regression for each SNP
* Adjust for covariates

### 5. Multiple Testing Correction

* Bonferroni or FDR

### 6. Visualization

* Manhattan plot
* QQ plot

## 2. Python Example Using PLINK and Pandas

We can use pandas, numpy, statsmodels, matplotlib for the analysis, and plink for genotype preprocessing.

* ### Step 1: Load Phenotype and Covariates

```
python
import pandas as pd

# Load phenotype file
pheno = pd.read_csv("phenotype.csv")  # columns: FID, IID, Trait
# Load covariates
covariates = pd.read_csv("covariates.csv")  # columns: FID, IID, Age, Sex, PC1, PC2
# Merge
data = pd.merge(pheno, covariates, on=['FID', 'IID'])
data.head()
```
* ### Step 2: Run QC in PLINK
```
# Filter SNPs with missing rate >5%, individuals with missing >5%
plink --bfile genotype --geno 0.05 --mind 0.05 --make-bed --out genotype_qc

# Filter SNPs by MAF >1%
plink --bfile genotype_qc --maf 0.01 --make-bed --out genotype_qc_maf

# Hardy-Weinberg filter (p<1e-6)
plink --bfile genotype_qc_maf --hwe 1e-6 --make-bed --out genotype_final
```
* ### Step 3: PCA for Population Stratification
  ```
  plink --bfile genotype_final --pca 10 --out pca
  ```
  ```
  pca = pd.read_csv("pca.eigenvec", delim_whitespace=True, header=None)
  pca.columns = ['FID', 'IID'] + [f'PC{i}' for i in range(1, 11)]
  data = pd.merge(data, pca, on=['FID', 'IID'])
  ```
* ### Step 4: GWAS Regression (Linear Regression Example)
  ```
  import numpy as np
  import statsmodels.api as sm

  # Load genotype dosage (e.g., from PLINK --recode A)
  geno = pd.read_csv("genotype_final.raw", delim_whitespace=True)  # columns: FID, IID, SNPs
  geno = geno.drop(['FID','IID'], axis=1)

  results = []
  for snp in geno.columns:
    X = pd.concat([geno[snp], data[['Age','Sex','PC1','PC2']]], axis=1)
    X = sm.add_constant(X)
    y = data['Trait']
    model = sm.OLS(y, X).fit()
    results.append([snp, model.params[snp], model.pvalues[snp]])

  gwas_df = pd.DataFrame(results, columns=['SNP','Beta','P'])
  gwas_df.head()
  ```
* ### Step 5: Multiple Testing Correction
  ```
  from statsmodels.stats.multitest import multipletests

  gwas_df['P_FDR'] = multipletests(gwas_df['P'], method='fdr_bh')[1]

  ```
* ### Step 6: Manhattan and QQ Plots
  ```
  import matplotlib.pyplot as plt
  import numpy as np

  # Manhattan plot
  plt.scatter(range(len(gwas_df)), -np.log10(gwas_df['P']), c='blue', s=2)
  plt.axhline(-np.log10(5e-8), color='red', linestyle='--')  # genome-wide significance
  plt.xlabel('SNPs')
  plt.ylabel('-log10(P)')
  plt.show()
  
  # QQ plot
  import scipy.stats as stats
  expected = -np.log10(np.linspace(1/len(gwas_df), 1, len(gwas_df)))
  observed = -np.log10(np.sort(gwas_df['P']))
  plt.scatter(expected, observed, s=5)
  plt.plot([0,max(expected)], [0,max(expected)], color='red', linestyle='--')
  plt.xlabel('Expected -log10(P)')
  plt.ylabel('Observed -log10(P)')
  plt.show()
  ```
  * ##### Notes

    * For large datasets, use specialized Python libraries like ```Hail``` or ```pandas-plink``` for efficiency.
    * Binary traits require logistic regression: ```sm.Logit(y, X).fit()```.
    * For covariates, always include PCs from PCA to correct for population stratification.
    * PLINK is still the most widely used tool for genotype QC and management.
