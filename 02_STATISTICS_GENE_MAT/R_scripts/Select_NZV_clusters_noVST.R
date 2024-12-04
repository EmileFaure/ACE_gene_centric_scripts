require(tidyverse) 
require(vegan)
require(caret)
#require(ggplot2)
#require(foreach)
#require(doParallel)

# Import Data
Gene_mat = readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T30GQ_transposed.rds")
MetaData=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/Metadata_GeneMat_KNN_MertzTF_1099Corr.rds")

# Select only CTD samples
#Gene_mat=Gene_mat[which(as.factor(row.names(Gene_mat)) %in% MetaData$ACE_seq_name),]
#rm(metadata_GeneMat_KNN)

# Remove nearzerovar
#ZV <- caret::nearZeroVar(Gene_mat, saveMetrics = TRUE, uniqueCut = 20)
#Gene_mat <- Gene_mat[,ZV$nzv == FALSE]
#rm(ZV)

#write.table(Gene_mat, file= "/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T30GQ_transposed_NZVuniquecut20.tsv",quote=FALSE,sep="\t")
#saveRDS(Gene_mat, "/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T30GQ_transposed_NZVuniquecut20.rds")

# Do the same for each separate size fraction
#smallsize = MetaData[which(MetaData$Size_fraction!=">3 µm"),1]
bigsize = MetaData[which(MetaData$Size_fraction==">3 µm"),1]

#Gene_mat_small=Gene_mat[which(as.factor(row.names(Gene_mat)) %in% smallsize),]
#ZV <- caret::nearZeroVar(Gene_mat_small, saveMetrics = TRUE, uniqueCut = 20)
#Gene_mat_small <- Gene_mat_small[,ZV$nzv == FALSE]
#rm(ZV)

#write.table(Gene_mat_small, file= "/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T30GQ_transposed_NZVuniquecut20_FreeLiving.tsv",quote=FALSE,sep="\t")
#saveRDS(Gene_mat_small, "/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T30GQ_transposed_NZVuniquecut20_FreeLiving.rds")
#rm(Gene_mat_small)

Gene_mat_big=Gene_mat[which(as.factor(row.names(Gene_mat)) %in% bigsize),]
ZV <- caret::nearZeroVar(Gene_mat_big, saveMetrics = TRUE, uniqueCut = 20)
Gene_mat_big <- Gene_mat_big[,ZV$nzv == FALSE]
rm(ZV)

write.table(Gene_mat_big, file= "/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T30GQ_transposed_NZVuniquecut20_Attached.tsv",quote=FALSE,sep="\t")
saveRDS(Gene_mat, "/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T30GQ_transposed_NZVuniquecut20_Attached.rds")
