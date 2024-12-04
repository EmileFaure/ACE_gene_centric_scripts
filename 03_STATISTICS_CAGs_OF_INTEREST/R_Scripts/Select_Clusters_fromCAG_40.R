library(tidyverse) 
library(data.table)

#CAG29
Gene_Mat=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T60MAXGQ_transposed_NZVuniquecut20_FL_RSquared10.rds")
AGC_list=fread("/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/CAG_40_FLRsquared10.txt",data.table=F,header=F)

# Select only AGC from the CAG
Gene_Mat <- Gene_Mat[,which(names(Gene_Mat) %in% AGC_list[,2])]

fwrite(Gene_Mat, file= "/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/GM_AGN_T60MAXGQ_transposed_CAG_40_FL.tsv",quote=FALSE,sep="\t",row.names=T)

