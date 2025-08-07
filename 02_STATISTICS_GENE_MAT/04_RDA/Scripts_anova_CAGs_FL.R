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
#The sign of some latent variables needs to be inversed to be interpreted as in the Landwehr et al. study
#Mail from Landwehr : "We have chosen to change the signs of some of the LVs that where produced by the algorithm, in order to make them align with certain physical processes"
#Also note that here LV are from 0 to 13, like in the data provided by Landwehr et al., need to add one to each LV number to match their study.
Envi_Mat$LV3_mean = -Envi_Mat$LV3_mean
Envi_Mat$LV12_mean = -Envi_Mat$LV12_mean
Envi_Mat$LV4_mean = -Envi_Mat$LV4_mean
Envi_Mat$LV11_mean = -Envi_Mat$LV11_mean
Envi_Mat$LV10_mean = -Envi_Mat$LV10_mean
Envi_Mat$LV5_mean = -Envi_Mat$LV5_mean
Envi_Mat$LV9_mean = -Envi_Mat$LV9_mean

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

