#Packages
library(data.table)
library(vegan)
library(tidyverse)
source("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/R_scripts/RDA/scoresRDA.R")

#Importing Data
#>10% R2
Gene_Mat=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_AGN_T60MAXGQ_transposed_NZVuniquecut20_FL_RSquared10.rds")
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

# PCA
PCA=rda(Gene_Mat_hel)
summary(PCA)

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

RDA=rda(Gene_Mat_hel ~ ., data = Envi_Mat_sel)
summary(RDA)
RsquareAdj(RDA) #77.5%

# Retrieve sample scores
scores=scoresRDA(RDA)

samples = scores$Sites
samples = as.data.frame(samples)
samples = merge(samples, Envi_Mat_sel, by = "row.names")
rownames(samples) = samples[,1]
samples = samples[,-1]
samples = merge(samples,Group_Var,by="row.names")
rownames(samples) = samples[,1]
samples = samples[,-1]

# Retrieve AGC scores
agcscore <- scores$Species
agcscore = as.data.frame(agcscore)

#Retrieve environmental scores
enviscore=scores$Biplot

#Plot the triplot
ggplot() +   geom_hline(yintercept = 0, linetype='dotted') +
  geom_vline(xintercept = 0, linetype='dotted') +
  labs(x = paste0("RDA1 (",
                      round(100 * scores$Eigenval.rel[1], 1), "%)"),
           y = paste0("RDA2 (",
                      round(100 * scores$Eigenval.rel[2], 1), "%)")) +
  theme(plot.title=element_text(hjust=0.5)) +
  geom_point(aes(x = agcscore$RDA1, y = agcscore$RDA2),size = 0.5, alpha = 0.4, col = "grey") +
  geom_point(aes(samples[,1], samples[,2], col = samples$Water_mass_simplified.x,shape=samples$MertzGlacier.x), size = 3, alpha = 0.9) +
  geom_segment(data=as.data.frame(enviscore),aes(xend = RDA1, yend = RDA2),x=0,y=0,size = 0.5, linetype="F1",color = 'steelblue4',arrow = arrow(length = unit(0.2,"cm"))) +
  geom_text(aes(enviscore[,1]*1.2, enviscore[,2]*1.2, label=rownames(enviscore)), color="steelblue4") +
  scale_colour_manual(values = levels(samples$Water_mass_simplified.col)) +
  labs(col = "Water mass")
ggsave("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/01-FIGURES/RDA_Triplot_FreeLiving_R2supp10.png",width=29.7,height=21,unit="cm")

# Select most informative environmental drivers
set.seed(1994)
rda.step.both = ordistep(rda(Gene_Mat_hel ~1,data=Envi_Mat_sel),scope=formula(RDA),direction="both",steps=100)
res.rda = rda(formula=formula(rda.step.both), data=Envi_Mat_sel)
saveRDS(res.rda,"/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/RDA/RDA_Object_PostVariableSelection_FreeLiving_R2supp10.rds")








