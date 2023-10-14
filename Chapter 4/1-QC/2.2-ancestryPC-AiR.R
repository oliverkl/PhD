### PCA with PC-AiR ###

# module load R

## load required packages
library(GENESIS)
library(SNPRelate)
library(SeqVarTools)
library(SeqArray)
library(Biobase)

## Create a GDS file from PLINK
## the plink file was merged with 1000g 
snpgdsBED2GDS( bed.fn = "${PLINKNAME}_1000gPCA_QC2.bed", fam.fn = "${PLINKNAME}_1000gPCA_QC2.fam", 
	bim.fn = "${PLINKNAME}_1000gPCA_QC2.bim", "${PLINKNAME}_1000gPCA_QC2.gds")

gdsfile <- "${PLINKNAME}_1000gPCA_QC2.gds"
gdsfmt::showfile.gds(closeall=TRUE) # make sure file is not already open
gds <- snpgdsOpen(gdsfile)

set.seed(100) # LD pruning has a random element; so make this reproducible
snpset <- snpgdsLDpruning(gds, method="corr", 
                          slide.max.bp=10e6, ld.threshold=sqrt(0.1))

pruned <- unlist(snpset, use.names=FALSE)

## King
king <- snpgdsIBDKING(gds, snp.id=pruned)

kingMat <- king$kinship
colnames(kingMat) <- rownames(kingMat) <- king$sample.id

## kinship mat
kinship <- snpgdsIBDSelection(king)

## PC-AiR
sampset <- pcairPartition(kinobj=kingMat, kin.thresh=2^(-9/2),
                          divobj=kingMat, div.thresh=-2^(-9/2))

sampset$unrels %>% data.frame() -> unrelated
names(unrelated) = "id"

# run PCA on unrelated set
pca.unrel <- snpgdsPCA(gds, sample.id=sampset$unrels, snp.id=pruned,eigen.cnt = 10)

# project values for relatives
snp.load <- snpgdsPCASNPLoading(pca.unrel, gdsobj=gds)

samp.load <- snpgdsPCASampLoading(snp.load, gdsobj=gds, sample.id=sampset$rels)

# combine unrelated and related PCs and order as in GDS file
pcs <- rbind(pca.unrel$eigenvect, samp.load$eigenvect)
rownames(pcs) <- c(pca.unrel$sample.id, samp.load$sample.id)

pc.df <- as.data.frame(pcs)
names(pc.df) <- 1:ncol(pcs)
pc.df$sample.id <- row.names(pcs)

saveRDS(pc.df, "pca.df.rds")