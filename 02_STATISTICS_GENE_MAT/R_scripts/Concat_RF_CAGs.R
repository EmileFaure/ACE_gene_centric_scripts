library(data.table)
library(tidyverse)

Gene_Mat=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T60MAXGQ_transposed_NZVuniquecut20_FL_RSquared10.rds")
RF_out=fread("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/RF_AGC_NZVuniquecut20_T60MAX_FL_All.txt",data.table=F)
CAGs=fread("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/CAGs_T60MAX_FL_RSquared10.txt",data.table=F,header=F)
names(CAGs)=c("CAG_ID","AGC_ID")
names(RF_out)[1]="AGC_ID"

# Gene Mat
Gene_Mat=as.data.frame(t(Gene_Mat))
Gene_Mat=merge(CAGs,Gene_Mat,by.x="AGC_ID",by.y="row.names")
Gene_Mat=Gene_Mat[,-1]
Gene_Mat=aggregate(.~CAG_ID,Gene_Mat,FUN=sum)
Gene_Mat=as.data.frame(t(Gene_Mat))
names(Gene_Mat)=Gene_Mat[1,]
Gene_Mat=Gene_Mat[-1,]
fwrite(Gene_Mat,"/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_CAGs_T60MAXGQ_transposed_NZVuniquecut20_FL_RSquared10.tsv",quote=FALSE,sep="\t",row.names=T)
rm(Gene_Mat)

# Results RF
RF_CAGs=merge(CAGs,RF_out,by="AGC_ID")
RF_CAGs=select(RF_CAGs,-c("AGC_ID"))
RF_CAGs_agg=aggregate(.~CAG_ID,RF_CAGs,FUN=mean)
RF_CAGs_SD=aggregate(R2_caret~CAG_ID,RF_CAGs,FUN=sd)
RF_CAGs_size=aggregate(R2_caret~CAG_ID,RF_CAGs,FUN=length)
RF_CAGs_agg$R2_caret_SD=RF_CAGs_SD$R2_caret
RF_CAGs_agg$CAG_Size=RF_CAGs_size$R2_caret
fwrite(RF_CAGs_agg,"/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/CAGs_T60MAX_FL_RSquared10_AggregateWithRF.txt",quote=FALSE,sep="\t")
