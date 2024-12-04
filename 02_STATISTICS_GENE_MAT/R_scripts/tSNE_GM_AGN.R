require(tidyverse) 
require(vegan)
require(caret)
require(Rtsne)
require(ggplot2)

# Import Data
#GroupingVar=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/metadata_metaG_grouping_variables.rds")
#metadata_GeneMat_KNN=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/metadata_GeneMat_KNN.rds")
Gene_Mat_AGN=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T30GQ_VST_transposed_NZVuniquecut20.rds")

# t-SNE
AGN_tSNE = Rtsne(Gene_Mat_AGN)
saveRDS(AGN_tSNE, "/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/tSNE_AGC_VST_NZVuniquecut20.rds")
