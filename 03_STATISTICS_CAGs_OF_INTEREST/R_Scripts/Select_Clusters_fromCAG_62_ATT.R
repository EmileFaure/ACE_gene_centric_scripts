library(tidyverse) 
library(data.table)

#CAG29
Gene_Mat=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T60MAXGQ_transposed_NZVuniquecut20_ATT_RSquared15.rds")
AGC_list=fread("/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/CAG_62_ATTRsquared15.txt",data.table=F,header=F)

# Select only AGC from the CAG
Gene_Mat <- Gene_Mat[,which(names(Gene_Mat) %in% AGC_list[,2])]

fwrite(Gene_Mat, file= "/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/GM_AGN_T60MAXGQ_transposed_CAG_62_ATT.tsv",quote=FALSE,sep="\t",row.names=T)

