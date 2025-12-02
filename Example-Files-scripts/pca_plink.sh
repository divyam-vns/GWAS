#!/bin/bash
# PCA for population stratification
plink --bfile ../results/genotype_final --pca 10 --out ../results/pca
