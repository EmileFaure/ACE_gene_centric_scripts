require(tidyverse) 
require(vegan)
require(caret)
require(ggplot2)
require(foreach)
require(doParallel)

# Import Data
Gene_mat = readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T30GQ_VST_transposed.rds")
metadata_GeneMat_KNN=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/metadata_GeneMat_KNN.rds")

# Select only CTD samples
Gene_mat=Gene_mat[which(as.factor(row.names(Gene_mat)) %in% metadata_GeneMat_KNN$ACE_seq_name),]
rm(metadata_GeneMat_KNN)

# Remove nearzerovar
ZV <- caret::nearZeroVar(Gene_mat, saveMetrics = TRUE, uniqueCut = 20)
Gene_mat <- Gene_mat[,ZV$nzv == FALSE]
rm(ZV)

write.table(Gene_mat, file= "/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T30GQ_VST_transposed_NZVuniquecut20.tsv",quote=FALSE,sep="\t")
saveRDS(Gene_mat, "/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T30GQ_VST_transposed_NZVuniquecut20.rds")

