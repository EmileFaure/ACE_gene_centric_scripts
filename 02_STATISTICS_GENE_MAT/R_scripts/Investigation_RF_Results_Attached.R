library(tidyverse)
library(data.table)
library(reshape2)

RF = fread("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/RF_AGC_NZVuniquecut20_T60MAX_ATT_All.txt",data.table=F)

# Visual exploration of results
ggplot() + 
  geom_density(aes(x=  RF$R2_caret), fill="deepskyblue3", col = "deepskyblue4", alpha = 0.6) +
  geom_segment(aes(x=mean(RF$R2_caret,na.rm=T), xend=mean(RF$R2_caret,na.rm=T), y=0, yend=2.6), linetype = "dashed", color = "darkblue") +
  labs(y = "Density", x = expression(R^{2} ~ "of cross-validated random forest models")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))

ggsave("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/01-FIGURES/Density_R2_Attached_T60MAX.pdf")

rank.rf.imp = as.data.frame(t(apply(as.matrix(RF[,c(6:ncol(RF))]), 1, rank)))
rank.rf.imp = 50-rank.rf.imp
boxplot = reshape2::melt(as.matrix(rank.rf.imp))
boxplot$Var2=sub("Importance.","",boxplot$Var2)
bpall = ggplot() + geom_boxplot(aes(x = reorder(Var2,value,FUN=median), y = value), data = boxplot) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face="bold", size=14), axis.text.y = element_text(size = 14), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
  labs(x = "", y = "Rank of importance in random forest regressions")
bpall
ggsave("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/01-FIGURES/Importance_Boxplot_Attached_T60MAX.png",width=29.7,height=21,unit="cm")

# Permutations
RF_Perm_1994 = fread("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/RF_AGC_TMAX60GQ_NZVuniquecut20_ATT_Permuted1994_All.txt",data.table=F)
RF_Perm_1997 = fread("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/RF_AGC_TMAX60GQ_NZVuniquecut20_ATT_Permuted1997_All.txt",data.table=F)

ggplot() +
  geom_density(aes(x=  RF$R2_caret), fill="deepskyblue3", col = "deepskyblue4", alpha = 0.6) +
  geom_segment(aes(x=mean(RF$R2_caret,na.rm=T), xend=mean(RF$R2_caret,na.rm=T), y=0, yend=24), linetype = "dashed", col = "deepskyblue4") +
  geom_density(aes(x=  RF_Perm_1994$R2_caret), fill="indianred1", col = "indianred2", alpha = 0.6) +  
  geom_segment(aes(x=quantile(RF_Perm_1994$R2_caret,na.rm=T,0.95), xend=quantile(RF_Perm_1994$R2_caret,na.rm=T,0.95), y=0, yend=24), linetype = "dashed", col = "indianred2") +
  geom_density(aes(x=  RF_Perm_1997$R2_caret), fill="indianred3", col = "indianred4", alpha = 0.6) +  
  geom_segment(aes(x=quantile(RF_Perm_1997$R2_caret,na.rm=T,0.95), xend=quantile(RF_Perm_1997$R2_caret,na.rm=T,0.95), y=0, yend=24), linetype = "dashed", col = "indianred4") +
  labs(y = "Density", x = expression(R^{2} ~ "of cross-validated random forest models")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))

ggsave("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/01-FIGURES/Density_R2_T60MAX_Attached_WithPerm_95percentile.pdf")

# Adding agnostos annotation context
AGC_cat=fread("/home/datawork-lmee-intranet-nos/ACE/04-GENE-CENTRIC/03-ANNOTATION/ACE_gene_cluster_files/cluster_repres_size_categories_AGC_HEADERS.tsv",data.table=F)

AGC_cat=AGC_cat[,-c(2,3)]

AGC_cat$cl_name=as.factor(AGC_cat$cl_name)
AGC_cat$category=as.factor(AGC_cat$category)
AGC_cat$category_singleton=as.factor(AGC_cat$category_singleton)

RF_AGC=merge(RF,AGC_cat,by.x="V1",by.y="cl_name",all.x=T,all.y=F)

RF_AGC$category=factor(RF_AGC$category,levels=c("SINGL","EU","GU","KWP","K"))
ggplot() +
  geom_density(aes(x=  RF_AGC$R2_caret, col= RF_AGC$category, fill=RF_AGC$category), alpha = 0.5) +
  labs(y = "Density", x = expression(R^{2} ~ "of cross-validated random forest models"), col = "AGNOSTOS cat.", fill = "AGNOSTOS cat.") +
  scale_colour_manual(values=c("hotpink2","#FF6666","#3399CC","#8fbc8f","#339966")) +
  scale_fill_manual(values=c("hotpink3","#FF6666","#3399CC","#8fbc8f","#339966")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))

ggsave("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/01-FIGURES/Density_R2_Attached_ByAGCCat_T60MAX.pdf")

