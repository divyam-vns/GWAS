#!/bin/bash
# QC pipeline using PLINK

# Step 1: Remove SNPs with missing rate >5% and individuals >5%
plink --bfile ../data/genotype --geno 0.05 --mind 0.05 --make-bed --out ../results/genotype_qc

# Step 2: Filter SNPs by minor allele frequency >1%
plink --bfile ../results/genotype_qc --maf 0.01 --make-bed --out ../results/genotype_qc_maf

# Step 3: Hardy-Weinberg equilibrium filter (p<1e-6)
plink --bfile ../results/genotype_qc_maf --hwe 1e-6 --make-bed --out ../results/genotype_final
