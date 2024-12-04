require(tidyverse) 
require(vegan)
require(caret)

# Import Data
Gene_mat=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T60MAXGQ_transposed.rds")
MetaData=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/Metadata_GeneMat_KNN_MertzTF_1099Corr.rds")

# Select only CTD samples
Gene_mat=Gene_mat[which(as.factor(row.names(Gene_mat)) %in% MetaData$ACE_seq_name),]

# Do the NZV selection for each separate size fraction
bigsize = MetaData[which(MetaData$Size_fraction==">3 µm"),1]

Gene_mat=Gene_mat[which(as.factor(row.names(Gene_mat)) %in% bigsize),]
ZV <- caret::nearZeroVar(Gene_mat, saveMetrics = TRUE, uniqueCut = 20)
Gene_mat <- Gene_mat[,ZV$nzv == FALSE]
rm(ZV)

saveRDS(Gene_mat, "/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_TMAX60GQ_transposed_NZVuniquecut20_Attached.rds")
write.table(Gene_mat, file= "/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_TMAX60GQ_transposed_NZVuniquecut20_Attached.tsv",quote=FALSE,sep="\t")
