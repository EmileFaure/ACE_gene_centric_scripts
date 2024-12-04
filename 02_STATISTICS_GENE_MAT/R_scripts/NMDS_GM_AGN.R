#require(tidyverse) 
require(vegan)

# Import Data
#Gene_Mat_AGN=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T30GQ_transposed_NZVuniquecut20.rds")
#Gene_Mat_AGN_small=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T30GQ_transposed_NZVuniquecut20_FreeLiving.rds")
Gene_Mat_AGN_big=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T30GQ_transposed_NZVuniquecut20_Attached.rds")

# NMDS
#AGN_NMDS = metaMDS(Gene_Mat_AGN, distance= "bray", trymax = 100, k = 2)
#saveRDS(AGN_NMDS, "/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/NMDS_AGC_NZVuniquecut20.rds")
#AGN_NMDS_small = metaMDS(Gene_Mat_AGN_small, distance= "bray", trymax = 100, k = 2)
#saveRDS(AGN_NMDS_small, "/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/NMDS_AGC_NZVuniquecut20_FreeLiving.rds")
AGN_NMDS_big = metaMDS(Gene_Mat_AGN_big, distance= "bray", trymax = 100, k = 2)
saveRDS(AGN_NMDS_big, "/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/NMDS_AGC_NZVuniquecut20_Attached.rds")
