#Packages
library(data.table)
library(vegan)
library(tidyverse)
library(data.table)
library(viridis)
source("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/R_scripts/RDA/scoresRDA.R")

#Importing Data
#>15% R2
Gene_Mat=fread("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_CAGs_T60MAXGQ_transposed_NZVuniquecut20_FL_RSquared10.tsv",data.table=F)
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

Group_Var=Group_Var[Group_Var$ACE_seq_name %in% row.names(Gene_Mat_hel),]
row.names(Group_Var)=Group_Var$ACE_seq_name
Group_Var=Group_Var[,-1]
#row match:
Group_Var=Group_Var[order(match(row.names(Group_Var),row.names(Gene_Mat_hel))),]

res.rda=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/RDA/RDA_Object_PostVariableSelection_CAGs_FreeLiving_R2supp10.rds")
RsquareAdj(res.rda) #72.9%
#Test significance (Warning, long to run):
set.seed(1994)
# Full model
anova.rda <- anova(res.rda) # Significant 0.001
anova.rda
# By axis
anova.rda.ax <- anova(res.rda, by = "axis", permutations = how(nperm = 499),parallel=8)
anova.rda.ax

