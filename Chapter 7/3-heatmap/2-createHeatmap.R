# R code to create heatmap for exploratory analysis

library(tidyverse)
library(ComplexHeatmap)

### Clustering heatmap

#1. Read in raw snp data
#2. Combine with phenotype data
#3. Heatmap

# read in clinical data
gabrg2_data <- readRDS("gefsplus_database.rds") %>% 
  filter(Gene == "GABRG2")

# read in raw genotype data from plink output
gaba.genotypes <- read_table("GABA_analysis/gabrg2_samps.gaba.snps.raw") %>% 
  select(-c(FID, SEX, MAT, PAT, PHENOTYPE))

# combine datasets
gaba.genotypes2 <- left_join(gabrg2_data, gaba.genotypes, by = "IID")

addmargins(table(gaba.genotypes2$severity, gaba.genotypes2$epilepsy, useNA = "ifany"))

# filter for p.R82Q carriers

gaba.genotypes6 <- gaba.genotypes5 %>% 
  filter(FID == "81650")

ped_ids <- read_tsv("GABRG2_fam/pedigree_IDs.txt")

gaba.genotypes7 <- left_join(gaba.genotypes6, ped_ids, by = c(IID = "DNA_N"))

# create matrix for heatmap

gaba_matrix <- gaba.genotypes7[1:32 , 8:36]

row.names <- gaba.genotypes7$Pedigree

epi.status <- as.data.frame(gaba.genotypes7$epilepsy)
row.names(epi.status) <- row.names
colnames(epi.status) <- "Epilepsy"
rownames(gaba_matrix) <- row.names

prs.group <- as.data.frame(gaba.genotypes7$allEpiPRSgroup)
row.names(prs.group) <- row.names
colnames(prs.group) <- "PRS"

gaba_matrix <- as.matrix(gaba_matrix)

### Heatmap ###

colors1 = structure(c("slategray1","dodgerblue1","blue"), names = c("0", "1", "2"))
colors2 = structure(c("firebrick2","grey57","green4"), names = c("high", "intermediate", "low"))
colors3 = structure(c("antiquewhite","mediumorchid4"), names = c("no", "yes"))


ht1 = Heatmap(gaba_matrix, name = "SNPs", col = colors1, column_title = "GABA receptor activity SNPs for GABRG2 p.Arg82Gln carriers",
              row_title = "Samples", rect_gp = gpar(col = "white", lwd = 1),
              column_names_rot = 60)
ht2 = Heatmap(prs.group, name = "PRS", col = colors2)
ht3 = Heatmap(epi.status, name = "Epilepsy", col = colors3)

ht1 + ht2 + ht3