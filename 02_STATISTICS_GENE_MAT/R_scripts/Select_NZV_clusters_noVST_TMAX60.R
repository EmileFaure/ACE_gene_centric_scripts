require(tidyverse) 
require(vegan)
require(caret)
require(data.table)
#require(ggplot2)
#require(foreach)
#require(doParallel)

# Import Genetic Data
#Gene_mat=fread("/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES/R_OBJECTS/GQ-Clust_TMAX_60PERC-Detec/Gene_Matrix_AGNOSTOS_T60MAXGQ_DESeq_Normalized_MCov.tsv")
#Corres_Samples=read.table("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/ACEsamples_With_Ace_seq_name.tsv",header=T,stringsAsFactors=T)

# Format gene matrix
#Gene_mat <- Gene_mat %>%
#    rename_at(vars(as.character(Corres_Samples$sample)), ~ as.character(Corres_Samples$ace_seq_name))

#Gene_mat <- as.data.frame(t(Gene_mat))
#names(Gene_mat)=Gene_mat[1,]
#Gene_mat=Gene_mat[-1,]

#Gene_mat[]<-data.frame(lapply(Gene_mat,as.numeric))

#saveRDS(Gene_mat, "/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T60MAXGQ_transposed.rds")
Gene_mat=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T60MAXGQ_transposed.rds")

# Import Envi Data
MetaData=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/Metadata_GeneMat_KNN_MertzTF_1099Corr.rds")

# Select only CTD samples
Gene_mat=Gene_mat[which(as.factor(row.names(Gene_mat)) %in% MetaData$ACE_seq_name),]
rm(MetaData)

# Remove nearzerovar
ZV <- caret::nearZeroVar(Gene_mat, saveMetrics = TRUE, uniqueCut = 20)
Gene_mat <- Gene_mat[,ZV$nzv == FALSE]
rm(ZV)

saveRDS(Gene_mat, "/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T60MAXGQ_transposed_NZVuniquecut20.rds")
write.table(Gene_mat, file= "/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T60MAXGQ_transposed_NZVuniquecut20.tsv",quote=FALSE,sep="\t")

