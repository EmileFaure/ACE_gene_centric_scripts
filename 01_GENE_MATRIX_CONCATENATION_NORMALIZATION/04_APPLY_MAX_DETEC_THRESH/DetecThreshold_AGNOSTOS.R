library(tidyverse)
library(data.table)

cov.mat <- fread("/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES/AGNOSTOS_CLSTRLVL_GENE_MAT_COV.tsv", sep="\t", header=TRUE,stringsAsFactors=TRUE)
det.mat <- fread("/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES/AGNOSTOS_CLSTRLVL_GENE_MAT_DET_Maximum.tsv", sep="\t", header=TRUE,stringsAsFactors=TRUE)

threshold=0.5

for(i in c(2:219)) {
 cov.mat[which(det.mat[[i]]<threshold),i]=0.0
}

write.table(cov.mat,file="/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES/AGNOSTOS_CLSTRLVL_GENE_MAT_COV_MAXIMUMDETECTHRESH_50percent.tsv",quote=F, sep="\t",row.names=F,dec=".",col.names=T)

threshold=0.6

for(i in c(2:219)) {
 cov.mat[which(det.mat[[i]]<threshold),i]=0.0
}

write.table(cov.mat,file="/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES/AGNOSTOS_CLSTRLVL_GENE_MAT_COV_MAXIMUMDETECTHRESH_60percent.tsv",quote=F, sep="\t",row.names=F,dec=".",col.names=T)

threshold=0.7

for(i in c(2:219)) {
 cov.mat[which(det.mat[[i]]<threshold),i]=0.0
}

write.table(cov.mat,file="/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES/AGNOSTOS_CLSTRLVL_GENE_MAT_COV_MAXIMUMDETECTHRESH_70percent.tsv",quote=F, sep="\t",row.names=F,dec=".",col.names=T)

threshold=0.8

for(i in c(2:219)) {
 cov.mat[which(det.mat[[i]]<threshold),i]=0.0
}

write.table(cov.mat,file="/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES/AGNOSTOS_CLSTRLVL_GENE_MAT_COV_MAXIMUMDETECTHRESH_80percent.tsv",quote=F, sep="\t",row.names=F,dec=".",col.names=T)

threshold=0.9

for(i in c(2:219)) {
 cov.mat[which(det.mat[[i]]<threshold),i]=0.0
}

write.table(cov.mat,file="/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES/AGNOSTOS_CLSTRLVL_GENE_MAT_COV_MAXIMUMDETECTHRESH_90percent.tsv",quote=F, sep="\t",row.names=F,dec=".",col.names=T)

threshold=1.0

for(i in c(2:219)) {
 cov.mat[which(det.mat[[i]]<threshold),i]=0.0
}

write.table(cov.mat,file="/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES/AGNOSTOS_CLSTRLVL_GENE_MAT_COV_MAXIMUMDETECTHRESH_100percent.tsv",quote=F, sep="\t",row.names=F,dec=".",col.names=T)
