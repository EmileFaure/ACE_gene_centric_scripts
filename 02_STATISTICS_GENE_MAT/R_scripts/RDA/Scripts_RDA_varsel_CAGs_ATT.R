#Packages
library(data.table)
library(vegan)
library(tidyverse)
library(data.table)
source("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/R_scripts/RDA/scoresRDA.R")

#Importing Data
#>15% R2
Gene_Mat=fread("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_CAGs_T60MAXGQ_transposed_NZVuniquecut20_ATT_RSquared15.tsv",data.table=F)
row.names(Gene_Mat)=Gene_Mat[,1]
Gene_Mat=Gene_Mat[,-1]
Envi_Mat=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/Metadata_GeneMat_KNN_MertzTF_1099Corr.rds")
Group_Var=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/metadata_metaG_grouping_variables.rds")

# Pre-investigations and processing

# Hellinger Transformation
Gene_Mat_hel=decostand(Gene_Mat,method="hellinger")

#Envi Data are already pre-processed (centered, scaled, knn imputed), just need to select and order rows
Envi_Mat_sel=Envi_Mat[Envi_Mat$ACE_seq_name %in% row.names(Gene_Mat_hel),]
row.names(Envi_Mat_sel)=Envi_Mat_sel$ACE_seq_name
Envi_Mat_sel=Envi_Mat_sel[,-1]
#row match:
Envi_Mat_sel=Envi_Mat_sel[order(match(row.names(Envi_Mat_sel),row.names(Gene_Mat_hel))),]

#Full RDA
RDA=rda(Gene_Mat_hel ~ ., data = Envi_Mat_sel)

# Select most informative environmental drivers
set.seed(1994)
rda.step.both = ordistep(rda(Gene_Mat_hel ~1,data=Envi_Mat_sel),scope=formula(RDA),direction="both",steps=100)
res.rda = rda(formula=formula(rda.step.both), data=Envi_Mat_sel)
saveRDS(res.rda,"/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/RDA/RDA_Object_PostVariableSelection_CAGs_Attached_R2supp15.rds")


