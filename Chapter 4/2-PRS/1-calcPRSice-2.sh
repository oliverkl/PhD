#!/bin/bash

#SBATCH --job-name=run_PRSice
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=3G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=oliver@wehi.edu.au

export APP=" "
export SUMSTATS=" "
export TARGET=" " #plink format files
export IN1=" "
export OUT1=" "
export PHENO=" "

module load R

PRSice()
{
Rscript PRSice.R \
--prsice ${APP}/PRSice_linux \
--base ${SUMSTATS} \
--A1 Allele1 \
--A2 Allele2 \
--stat Beta \
--beta \
--binary-target T \
--bp BP \
--chr CHR \
--pvalue P-value \
--snp MarkerName \
--target ${IN1}/${TARGET} \
--out ${OUT1}/${PHENO} \
--bar-levels 0.1,0.5 \
--fastscore \
--no-regress \
--no-full \
--model add \
--missing MEAN_IMPUTE \
--score avg \
--clump-kb 250kb \
--clump-p 1.0 \
--clump-r2 0.1 \
--upper 1 \
--print-snp
}

PRSice
