setwd("path/to/working/directory")
rm(list = ls())

library(openxlsx)
library(data.table)
library(org.Hs.eg.db)
# library(dplyr)
# library(tidyr)
library(clusterProfiler)
# source('fun_CalculateHartEssentialGenes_gmcsTH.R')

th <- 0.05


CCLE_Exp <- as.data.frame(fread("DepMap_23Q4/OmicsExpressionProteinCodingGenesTPMLogp1.csv"))
colnames(CCLE_Exp) <- gsub("[\\(\\)]", "", regmatches(colnames(CCLE_Exp), gregexpr("\\(.*?\\)", colnames(CCLE_Exp))))
genes <- colnames(CCLE_Exp)
CCLE_Exp <- as.data.frame(t(CCLE_Exp))

gc()
colnames(CCLE_Exp) <- CCLE_Exp[1,]
CCLE_Exp$ENTREZID <- genes
toENSEMBL <- clusterProfiler::bitr(CCLE_Exp[,"ENTREZID"], "ENTREZID", "ENSEMBL", org.Hs.eg.db)
CCLE_Exp <- merge(toENSEMBL, CCLE_Exp, by = "ENTREZID")
CCLE_Exp <- CCLE_Exp[!duplicated(CCLE_Exp$ENTREZID),]
CCLE_Exp <- CCLE_Exp[!duplicated(CCLE_Exp$ENSEMBL),]
rownames(CCLE_Exp) <- CCLE_Exp[,2]
CCLE_Exp <- CCLE_Exp[,-c(1,2)]

CCLE_genes <- rownames(CCLE_Exp)

CCLE_Exp<- apply(CCLE_Exp,2,as.numeric)
rownames(CCLE_Exp) <- CCLE_genes

gMCSs_Human <- read.xlsx("calculated_gMCSs_gMIS.xlsx", sheet = "Human1")


gMCSs.ENSEMBL <- as.matrix(gMCSs_Human)
rownames(gMCSs.ENSEMBL) <- as.character(seq(1,nrow(gMCSs.ENSEMBL)))
gMCSs.ENSEMBL <- unique(gMCSs.ENSEMBL)
gMCSs.ENSEMBL[is.na(gMCSs.ENSEMBL)] <- ""

genes_in <- apply(gMCSs.ENSEMBL, 1, function(x){
  x %in% CCLE_genes
})

genes_in <- colSums(genes_in) 

lengths <- apply(gMCSs.ENSEMBL, 2, function(x){
  pos <- which(x == "")
  pos[length(pos)]
})

lengths <- do.call(cbind,lengths)
lengths[length(lengths)+1] <- nrow(gMCSs.ENSEMBL)

a <- 0
keep <- array()
for (i in 1:length(lengths)){
  pos <- lengths[i]
  keep[(a+1):(pos+a)] <- genes_in[(a+1):(pos+a)] == i
  a <- pos
}

gMCSs.ENSEMBL <- gMCSs.ENSEMBL[which(keep),]
gMCSs.ENSEMBL.4.aux <- apply(gMCSs.ENSEMBL,1,function(x){paste0(x[x!=""],collapse = "--")})

gMCSs.ENSEMBL.list <- strsplit(gMCSs.ENSEMBL.4.aux, '--')

CCLE_Exp[is.na(CCLE_Exp)] <- 0 
gene.exp <- CCLE_Exp
dim(gene.exp)

exp_max <- matrix(0, nrow = length(gMCSs.ENSEMBL.list), ncol = ncol(CCLE_Exp))
names_max2 <- matrix("", nrow = length(gMCSs.ENSEMBL.list), ncol = ncol(CCLE_Exp))

CCLE_Exp[is.na(CCLE_Exp)] <- 0

for(i in seq_along(gMCSs.ENSEMBL.list)) {
  genes <- gMCSs.ENSEMBL.list[[i]]
  exp <- CCLE_Exp[genes, , drop = FALSE]
  if (length(genes) == 1) {
    exp_max[i, ] <- exp
    names_max2[i, ] <- genes
  } else {
    exp_max[i, ] <- apply(exp, 2, max)
    names_max <- apply(exp, 2, function(col) genes[which.max(col)])
    names_max2[i, ] <- names_max
  }
}

gc()

threshold <- array()
for (i in 1:ncol(names_max2)){
  unique_genes <- !duplicated(names_max2[,i])
  aux <- exp_max[unique_genes,i]
  threshold[i] <- quantile(aux[aux>0], th)
}

threshold <- as.data.frame(threshold)
rownames(threshold) <- colnames(CCLE_Exp)

write.table(threshold, "th_log2TPM_CCLE.txt")

gc()
