#!/bin/bash

# Harmonise datasets to HRC reference ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz 
# using GenotypeHarmonizer: https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer

# set environmental variables

# directories
export IN1=" "
export OUT1=" "
export REF=" "

# file name
export PLINK1=" "


java -jar GenotypeHarmonizer-1.4.20-SNAPSHOT/GenotypeHarmonizer.jar \
--input ${OUT1}/${PLINK1}_QC4_Euro \
--inputType PLINK_BED \
--ref ${REF}/HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz \
--refType VCF \
--update-id \
--update-reference-allele \
--output ${OUT1}/${PLINK1}_QC5_HRC-aligned \
--outputType PLINK_BED

## Prepare harmonised dataset for upload to MIS

MAKEINPUT()
{
  for chr in {1..23}
  do
  plink \
  --bfile ${OUT1}/${PLINK1}_QC5_HRC-aligned \
  --chr ${chr} \
  --chr-output M \
  --set-hh-missing \
  --allow-no-sex \
  --recode vcf bgz \
  --out ${OUT1}/MIS/${PLINK1}_QC5_HRC-aligned_chr${chr}
  done
}

MAKEINPUT