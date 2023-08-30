#!/bin/bash

# set environmental variables

# directories
export IN1=" "
export OUT1=" "

# file name
export PLINK1=" "

module load plink

# remove SNPs with call rate <98%
plink \
--bfile ${IN1}/${PLINK1} \
--allow-no-sex \
--geno 0.02 \
--make-bed \
--out ${OUT1}/${PLINK1}_QC1_snp_0.98

# remove samples with call rate <98%
plink \
--bfile ${OUT1}/${PLINK1}_QC1_snp_0.98 \
--allow-no-sex \
--mind 0.02 \
--make-bed \
--out ${OUT1}/${PLINK1}_QC2_samp_0.98

# remove SNPs with MAF <0.5% or HWE deviations
plink \
--allow-no-sex \
--bfile ${OUT1}/${PLINK1}_QC2_samp_0.98 \
--maf 0.05 \
--hwe 1e-6 midp \
--make-bed  \
--out ${OUT1}/${PLINK1}_QC3

# identify outliers in heterozygosity, sex failures and calculate pairwise IBD
plink \
--bfile ${OUT1}/${PLINK1}_QC3 \
--het \
--check-sex \
--check-sex ycount \
--genome
