require(tidyverse) 
require(vegan)
require(data.table)

# Import Data
#Gene_Mat_AGN_small=fread("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_TMAX60GQ_transposed_NZVuniquecut20_FreeLiving.tsv", data.table=F)
Gene_Mat_AGN_big=fread("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_TMAX60GQ_transposed_NZVuniquecut20_Attached.tsv", data.table=F)
names(Gene_Mat_AGN_big)[1]="ACE_seq_name"

MetaData = readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/Metadata_GeneMat_KNN_MertzTF_1099Corr.rds")
MetaData <- MetaData %>% select(ACE_seq_name,Water_mass_simplified)

#Gene_Mat_AGN_small <- merge(Gene_Mat_AGN_small,MetaData,by="ACE_seq_name")
Gene_Mat_AGN_big <- merge(Gene_Mat_AGN_big,MetaData,by="ACE_seq_name")

# ANOSIM
#ANOSIM_small = anosim(Gene_Mat_AGN_small[,-c(1,ncol(Gene_Mat_AGN_small))], Gene_Mat_AGN_small[,"Water_mass_simplified"], distance = "bray", permutations = 9999)
#summary(ANOSIM_small)
#saveRDS(ANOSIM_small, "/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/ANOSIM_AGC_NZVuniquecut20_T60MAX_FreeLiving_WaterMasses.rds")

ANOSIM_big = anosim(Gene_Mat_AGN_big[,-c(1,ncol(Gene_Mat_AGN_big))], Gene_Mat_AGN_big[,"Water_mass_simplified"], distance = "bray", permutations = 9999)
summary(ANOSIM_big)
saveRDS(ANOSIM_big, "/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/ANOSIM_AGC_NZVuniquecut20_T60MAX_Attached_WaterMasses.rds")
