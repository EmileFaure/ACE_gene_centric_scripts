library(data.table)

RF = fread("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/RF_AGC_NZVuniquecut20_T60MAX_FL_All.txt",data.table=F)

Gene_Mat_AGN=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_TMAX60GQ_transposed_NZVuniquecut20_FreeLiving.rds")

RF_10=RF$V1[which(RF$R2_caret>0.1)]
#write.table(RF_10,"/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/AGC_CLUSTERS_RSQUARED10_FL.txt", quote=F, row.names=F)
RF_25=RF$V1[which(RF$R2_caret>0.25)]
#write.table(RF_25,"/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/AGC_CLUSTERS_RSQUARED25_FL.txt", quote=F, row.names=F)
RF_50=RF$V1[which(RF$R2_caret>0.5)]
#write.table(RF_50,"/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/AGC_CLUSTERS_RSQUARED50_FL.txt", quote=F, row.names=F)

Gene_Mat_AGN_sel=Gene_Mat_AGN[,which(names(Gene_Mat_AGN) %in% RF_10)]
#saveRDS(Gene_Mat_AGN_sel,"/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T60MAXGQ_transposed_NZVuniquecut20_FL_RSquared10.rds")
fwrite(Gene_Mat_AGN_sel,"/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T60MAXGQ_transposed_NZVuniquecut20_FL_RSquared10.tsv",quote=FALSE,sep="\t",row.names=TRUE)
Gene_Mat_AGN_sel=Gene_Mat_AGN[,which(names(Gene_Mat_AGN) %in% RF_25)]
fwrite(Gene_Mat_AGN_sel,"/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T60MAXGQ_transposed_NZVuniquecut20_FL_RSquared25.tsv",quote=FALSE,sep="\t",row.names=TRUE)
#saveRDS(Gene_Mat_AGN,"/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T60MAXGQ_transposed_NZVuniquecut20_FL_RSquared25.rds")
Gene_Mat_AGN_sel=Gene_Mat_AGN[,which(names(Gene_Mat_AGN) %in% RF_50)]
#saveRDS(Gene_Mat_AGN_sel,"/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T60MAXGQ_transposed_NZVuniquecut20_FL_RSquared50.rds")
fwrite(Gene_Mat_AGN_sel,"/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T60MAXGQ_transposed_NZVuniquecut20_FL_RSquared50.tsv",quote=FALSE,sep="\t",row.names=TRUE)
