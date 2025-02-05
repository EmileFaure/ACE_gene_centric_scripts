### Concatenation and normalization of Anvio-profile-blitz outputs ###

library(tidyverse)
library(DESeq2)
library(readxl)
library(data.table)
library(stringr)
#setwd("/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES")


###### Normalizing Mean Coverages

clust.matrix.cov <- fread("/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES/AGNOSTOS_GQ_CLSTRLVL_GENE_MAT_COV_MAXIMUMDETECTHRESH_60percent.tsv", sep="\t", header=TRUE,stringsAsFactors=TRUE)
clust.matrix.cov=as.data.frame(clust.matrix.cov)
#saveRDS(clust.matrix.cov,file="/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES/R_OBJECTS/Gene_Matrix_AGNOSTOS_MCov.RDS")

#Adding AGC_ in ID name
clust.matrix.cov[,1]=paste("AGC",clust.matrix.cov[,1],sep="_")

# Using DESeq2 approach

# Need to transform decimal values to integers
clust.matrix.cov <- clust.matrix.cov %>% mutate_if(is.numeric,round)
#saveRDS(clust.matrix.cov, file="/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/R_OBJECTS/07_GENE_MATRICES/Gene_Matrix_CDHit_ROUNDED.RDS")

coldata <- data.frame(samples = names(clust.matrix.cov)[-1])
row.names(coldata) = coldata$samples

dds <- DESeqDataSetFromMatrix(countData = clust.matrix.cov, colData = coldata, design= ~ 1, tidy=TRUE)
dds <- estimateSizeFactors(dds)
saveRDS(dds, file="/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES/R_OBJECTS/GQ-Clust_TMAX_60PERC-Detec/Gene_Matrix_AGNOSTOS_T60MAXGQ_DESeqObject.RDS")

#dds <- readRDS("/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES/R_OBJECTS/Gene_Matrix_AGNOSTOS_DESeqObject.RDS")

size_factors <- as.data.frame(sizeFactors(dds))
saveRDS(size_factors, file="/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES/R_OBJECTS/GQ-Clust_TMAX_60PERC-Detec/Size_Factors_DESeq_AGNOSTOS_T60MAXGQ.RDS")
write.table(size_factors,"/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES/R_OBJECTS/GQ-Clust_TMAX_60PERC-Detec/Size_Factors_DESeq_AGNOSTOS_T60MAXGQ.tsv",quote=FALSE,sep="\t")
rm(size_factors)

DS.normalized_clust.matrix <- counts(dds, normalized=TRUE)
saveRDS(DS.normalized_clust.matrix, file="/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES/R_OBJECTS/GQ-Clust_TMAX_60PERC-Detec/Gene_Matrix_AGNOSTOS_T60MAXGQ_DESeq_Normalized_MCov.RDS")
write.table(DS.normalized_clust.matrix,"/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES/R_OBJECTS/GQ-Clust_TMAX_60PERC-Detec/Gene_Matrix_AGNOSTOS_T60MAXGQ_DESeq_Normalized_MCov.tsv",quote=FALSE,sep="\t")
rm(DS.normalized_clust.matrix)

DS.vst_clust.matrix <- varianceStabilizingTransformation(dds, blind=TRUE, fitType="parametric")
DS.vst_clust.matrix <- as.data.frame(assay(DS.vst_clust.matrix))
saveRDS(DS.vst_clust.matrix, file="/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES/R_OBJECTS/GQ-Clust_TMAX_60PERC-Detec/Gene_Matrix_AGNOSTOS_T60MAXGQ_DESeq_Normalized_VST_MCov.RDS")
write.table(DS.vst_clust.matrix,"/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES/R_OBJECTS/GQ-Clust_TMAX_60PERC-Detec/Gene_Matrix_AGNOSTOS_T60MAXGQ_DESeq_Normalized_VST_MCov.tsv",quote=FALSE,sep="\t")
rm(dds,DS.vst_clust.matrix)



















