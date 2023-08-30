#!/bin/bash

### Harmonise to 1000 Genomes data prior to PCA ###

# set environmental variables
export DIRIN=" "
export DIRQC=" "
export PLINKNAME=" "
export REF=" "

java -jar GenotypeHarmonizer-1.4.20-SNAPSHOT/GenotypeHarmonizer.jar \
--input ${DIRIN}/${PLINKNAME} \
--inputType PLINK_BED \
--ref ${REF}/1000g_Hg19/1000g_ALLchr.phase3_shapeit2_g0.05_m0.01_h10-15_onlySNP \
--refType PLINK_BED \
--update-id \
--update-reference-allele \
--output ${DIRQC}/${PLINKNAME}_1000g_aligned \
--outputType PLINK_BED

# merge our dataset with 1000 Genomes
plink \
--bfile ${DIRQC}/${PLINKNAME}_1000g_aligned \
--extract range ${DIRQC}/extract-range_${PLINKNAME}_1000g_aligned-SNPs.txt \
--bmerge ${REF}/1000g_Hg19/1000g_ALLchr.phase3_shapeit2_g0.05_m0.01_h10-15_onlySNP \
--make-bed \
--out ${DIRQC}/${PLINKNAME}_1000gPCA_QC0

# standard QC on merged dataset
plink \
--bfile ${DIRQC}/${PLINKNAME}_1000gPCA_QC0 \
--chr 1-22 \
--geno 0.02 \
--maf 0.05 \
--mind 0.02 \
--hwe 1e-6 midp \
--make-bed \
--out ${DIRQC}/${PLINKNAME}_1000gPCA_QC1

## prune to set of independent SNPs
plink \
--allow-no-sex \
--bfile ${DIRQC}/${PLINKNAME}_1000gPCA_QC1 \
--indep-pairwise 2000 100 0.1 \
--out ${DIRQC}/${PLINKNAME}_1000gPCA_QC2

plink \
--bfile ${DIRQC}/${PLINKNAME}_1000gPCA_QC1 \
--extract ${DIRQC}/${PLINKNAME}_1000gPCA_QC2.prune.in \
--make-bed \
--out ${DIRQC}/${PLINKNAME}_1000gPCA_QC2
