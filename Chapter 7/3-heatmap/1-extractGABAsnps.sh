#!/bin/bash

#Plink to select GABRG2 samples and GABA SNPs

plink \
--bfile orig.data \
--keep gabrg2_samps.txt \
--out gabrg2_samps \
--allow-no-sex \
--make-bed 

# Filter SNPs to those with missingness <0.5%
plink \
--bfile gabrg2_samps \
--geno 0.005 \
--out gabrg2_samps_snp0.995 \
--allow-no-sex \
--make-bed 

#Extract all epilepsy PRS gaba SNPs and convert to raw format
plink \
--bfile gabrg2_samps_snp0.995 \
--extract gaba_snps.txt \
--out gabrg2_samps.gaba.snps \
--allow-no-sex \
--recodeA 
