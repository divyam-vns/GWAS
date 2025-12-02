#!/usr/bin/env python3
import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

# Load phenotype and covariates
pheno = pd.read_csv('../data/phenotype.csv')
cov = pd.read_csv('../data/covariates.csv')
data = pd.merge(pheno, cov, on=['FID','IID'])

# Load genotype (PLINK --recode A output)
geno = pd.read_csv('../results/genotype_final.raw', delim_whitespace=True)
geno = geno.drop(['FID','IID'], axis=1)

# Merge PCA if needed
pca = pd.read_csv('../results/pca.eigenvec', delim_whitespace=True, header=None)
pca.columns = ['FID','IID'] + [f'PC{i}' for i in range(1,11)]
data = pd.merge(data, pca, on=['FID','IID'])

# GWAS: Linear regression for each SNP
results = []
for snp in geno.columns:
    X = pd.concat([geno[snp], data[['Age','Sex','PC1','PC2']]], axis=1)
    X = sm.add_constant(X)
    y = data['Trait']
    model = sm.OLS(y, X).fit()
    results.append([snp, model.params[snp], model.pvalues[snp]])

gwas_df = pd.DataFrame(results, columns=['SNP','Beta','P'])
gwas_df['P_FDR'] = multipletests(gwas_df['P'], method='fdr_bh')[1]

# Save results
gwas_df.to_csv('../results/gwas_results.csv', index=False)

# Manhattan plot
plt.scatter(range(len(gwas_df)), -np.log10(gwas_df['P']), c='blue', s=2)
plt.axhline(-np.log10(5e-8), color='red', linestyle='--')
plt.xlabel('SNPs')
plt.ylabel('-log10(P)')
plt.title('Manhattan Plot')
plt.savefig('../results/manhattan_plot.png', dpi=300)
plt.close()

# QQ plot
expected = -np.log10(np.linspace(1/len(gwas_df),1,len(gwas_df)))
observed = -np.log10(np.sort(gwas_df['P']))
plt.scatter(expected, observed, s=5)
plt.plot([0,max(expected)], [0,max(expected)], color='red', linestyle='--')
plt.xlabel('Expected -log10(P)')
plt.ylabel('Observed -log10(P)')
plt.title('QQ Plot')
plt.savefig('../results/qq_plot.png', dpi=300)
plt.close()
