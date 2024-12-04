require(tidyverse) 
require(vegan)
require(data.table)

# Import Data
#Gene_Mat_AGN_small=fread("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_TMAX60GQ_transposed_NZVuniquecut20_FreeLiving.tsv", data.table=F)
#names(Gene_Mat_AGN_small)[1]="ACE_seq_name"
Gene_Mat_AGN_big=fread("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_TMAX60GQ_transposed_NZVuniquecut20_Attached.tsv", data.table=F)
names(Gene_Mat_AGN_big)[1]="ACE_seq_name"

MetaData = readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/Metadata_GeneMat_KNN_MertzTF_1099Corr.rds")
MetaData <- MetaData %>% select(ACE_seq_name,Water_mass_simplified)

#Gene_Mat_AGN_small <- merge(Gene_Mat_AGN_small,MetaData,by="ACE_seq_name")
Gene_Mat_AGN_big <- merge(Gene_Mat_AGN_big,MetaData,by="ACE_seq_name")

# ADONIS
#ADONIS_small = adonis(Gene_Mat_AGN_small[,-c(1,ncol(Gene_Mat_AGN_small))]~Water_mass_simplified,data=Gene_Mat_AGN_small, method = "bray", perm = 999)
#print(ADONIS_small)
#saveRDS(ADONIS_small, "/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/ADONIS_AGC_NZVuniquecut20_T60MAX_FreeLiving_WaterMasses.rds")

ADONIS_big = adonis(Gene_Mat_AGN_big[,-c(1,ncol(Gene_Mat_AGN_big))]~Water_mass_simplified,data=Gene_Mat_AGN_big, method = "bray", perm = 999)
print(ADONIS_big)
saveRDS(ADONIS_big, "/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/ADONIS_AGC_NZVuniquecut20_T60MAX_Attached_WaterMasses.rds")
