#Packages
library(data.table)
library(vegan)
library(tidyverse)
library(data.table)
library(viridis)
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

# PCA
#PCA=rda(Gene_Mat_hel)
#summary(PCA)

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

#RDA=rda(Gene_Mat_hel ~ ., data = Envi_Mat_sel)
#summary(RDA)
#RsquareAdj(RDA) #80.45%

# Retrieve sample scores
#scores=scoresRDA(RDA)

#samples = scores$Sites
#samples = as.data.frame(samples)
#samples = merge(samples, Envi_Mat_sel, by = "row.names")
#rownames(samples) = samples[,1]
#samples = samples[,-1]
#samples = merge(samples,Group_Var,by="row.names")
#rownames(samples) = samples[,1]
#samples = samples[,-1]

# Retrieve AGC scores
#agcscore <- scores$Species
#agcscore = as.data.frame(agcscore)

#enviscore=scores$Biplot

#Plot the triplot
#ggplot() +   geom_hline(yintercept = 0, linetype='dotted') +
#  geom_vline(xintercept = 0, linetype='dotted') +
#  labs(x = paste0("RDA1 (",
#                      round(100 * scores$Eigenval.rel[1], 1), "%)"),
#           y = paste0("RDA2 (",
#                      round(100 * scores$Eigenval.rel[2], 1), "%)")) +
#  theme(plot.title=element_text(hjust=0.5)) +
#  geom_point(aes(x = agcscore$RDA1, y = agcscore$RDA2),size = 0.5, alpha = 0.4, col = "grey") +
#  geom_point(aes(samples[,1], samples[,2], col = samples$Water_mass_simplified.x,shape=samples$MertzGlacier.x), size = 3, alpha = 0.9) +
#  geom_segment(data=as.data.frame(enviscore),aes(xend = RDA1, yend = RDA2),x=0,y=0,size = 0.5, linetype="F1",color = 'steelblue4',arrow = arrow(length = unit(0.2,"cm"))) +
#  geom_text(aes(enviscore[,1]*1.2, enviscore[,2]*1.2, label=rownames(enviscore)), color="steelblue4") +
#  scale_colour_manual(values = levels(samples$Water_mass_simplified.col)) +
#  labs(col = "Water mass") +
#  theme_bw()
#ggsave("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/01-FIGURES/RDA_Triplot_CAGs_Attached_R2supp15.png",width=29.7,height=21,unit="cm")

#Check R-square value for the clusters in the triplot
RF=fread("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/CAGs_T60MAX_ATT_RSquared15_AggregateWithRF.txt",data.table=F)
RF=select(RF,c("CAG_ID","R2_caret","R2_caret_SD","CAG_Size"))
#agcscore=merge(agcscore,RF,by.x="row.names",by.y="CAG_ID")
#ggplot() +   geom_hline(yintercept = 0, linetype='dotted') +
#  geom_vline(xintercept = 0, linetype='dotted') +
#  labs(x = paste0("RDA1 (",
#                      round(100 * scores$Eigenval.rel[1], 1), "%)"),
#           y = paste0("RDA2 (",
#                      round(100 * scores$Eigenval.rel[2], 1), "%)")) +
#  theme(plot.title=element_text(hjust=0.5)) +
#  geom_point(aes(x = agcscore$RDA1, y = agcscore$RDA2, col=agcscore$R2_caret,size=agcscore$CAG_Size), alpha = 0.6) +
#  scale_color_viridis(option="plasma") +
#  theme_bw()
#ggsave("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/01-FIGURES/RDA_CAGplot_ATT_R2supp15_colbyR2.png",width=21,height=30,unit="cm")
# The plot of only samples and envi to go with it :
#ggplot() +   geom_hline(yintercept = 0, linetype='dotted') +
#  geom_vline(xintercept = 0, linetype='dotted') +
#  labs(x = paste0("RDA1 (",
#                      round(100 * scores$Eigenval.rel[1], 1), "%)"),
#           y = paste0("RDA2 (",
#                      round(100 * scores$Eigenval.rel[2], 1), "%)")) +
#  theme(plot.title=element_text(hjust=0.5)) +
#  geom_segment(data=as.data.frame(enviscore),aes(xend = RDA1, yend = RDA2),x=0,y=0,size = 0.5, linetype="F1",color = 'steelblue4',arrow = arrow(length = unit(0.2,"cm"))) +
#  geom_text(aes(enviscore[,1]*1.2, enviscore[,2]*1.2, label=rownames(enviscore)), color="steelblue4") +
#  geom_point(aes(samples[,1], samples[,2], col = samples$Water_mass_simplified.x,shape=samples$MertzGlacier.x), size = 3, alpha = 0.9) +
#  scale_colour_manual(values = levels(samples$Water_mass_simplified.col)) +
#  labs(col = "Water mass") +
#  theme_bw()
#ggsave("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/01-FIGURES/RDA_CAGplot_ATT_R2supp15_onlyEnviSites.png",width=21,height=30,unit="cm")

# Select most informative environmental drivers : Run varsel script (very long)
#set.seed(1994)
#rda.step.both = ordistep(rda(Gene_Mat_hel ~1,data=Envi_Mat_sel),scope=formula(RDA),direction="both",steps=100)
#res.rda = rda(formula=formula(rda.step.both), data=Envi_Mat_sel)
#saveRDS(res.rda,"/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/RDA/RDA_Object_PostVariableSelection_CAGs_FreeLiving_R2supp10.rds")
res.rda=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/RDA/RDA_Object_PostVariableSelection_CAGs_Attached_R2supp15.rds")
RsquareAdj(res.rda) #74.6%
#Test significance (Warning, long to run):
#set.seed(1994)
# Full model
#anova.rda <- anova(res.rda) # Significant 0.001
# By axis
#anova.rda.ax <- anova(res.rda, by = "axis", permutations = how(nperm = 499)) # 11 significant axes

# Retrieve sample scores
scores=scoresRDA(res.rda)

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
  labs(col = "Water mass") +
  theme_bw()
ggsave("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/01-FIGURES/RDA_Triplot_CAGs_Attached_R2supp15_PostVarSel.pdf",width=23,height=23,unit="cm")

#Check R-square value for the clusters in the triplot
agcscore=merge(agcscore,RF,by.x="row.names",by.y="CAG_ID")
ggplot() +   geom_hline(yintercept = 0, linetype='dotted') +
  geom_vline(xintercept = 0, linetype='dotted') +
  labs(x = paste0("RDA1 (",
                      round(100 * scores$Eigenval.rel[1], 1), "%)"),
           y = paste0("RDA2 (",
                      round(100 * scores$Eigenval.rel[2], 1), "%)")) +
  theme(plot.title=element_text(hjust=0.5)) +
  geom_point(aes(x = agcscore$RDA1, y = agcscore$RDA2, col=agcscore$R2_caret,size=agcscore$CAG_Size), alpha = 0.6) +
  scale_color_viridis(option="plasma") +
  theme_bw()
ggsave("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/01-FIGURES/RDA_CAGplot_ATT_R2supp15_colbyR2_PostVarSel.pdf",width=23,height=23,unit="cm")

# The plot of only samples and envi to go with it :
ggplot() +   geom_hline(yintercept = 0, linetype='dotted') +
  geom_vline(xintercept = 0, linetype='dotted') +
  labs(x = paste0("RDA1 (",
                      round(100 * scores$Eigenval.rel[1], 1), "%)"),
           y = paste0("RDA2 (",
                      round(100 * scores$Eigenval.rel[2], 1), "%)")) +
  theme(plot.title=element_text(hjust=0.5)) +
  geom_segment(data=as.data.frame(enviscore),aes(xend = RDA1, yend = RDA2),x=0,y=0,size = 0.5, linetype="F1",color = 'steelblue4',arrow = arrow(length = unit(0.2,"cm"))) +
  geom_text(aes(enviscore[,1]*1.05, enviscore[,2]*1.05, label=rownames(enviscore)), color="steelblue4") +
  geom_point(aes(samples[,1], samples[,2], col = samples$Water_mass_simplified.x,shape=samples$MertzGlacier.x), size = 3, alpha = 0.9) +
  scale_colour_manual(values = levels(samples$Water_mass_simplified.col)) +
  labs(col = "Water mass") +
  theme_bw()
ggsave("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/01-FIGURES/RDA_CAGplot_ATT_R2supp15_onlyEnviSites_PostVarSel.pdf",width=23,height=23,unit="cm")

#Identification of key CAC
# The ones obviously linked with mertz
agcscore[which(agcscore$RDA1<(-0.5)),]
#CAG_15 -2.303489 -0.3638991 0.6107939  0.11681604   177373
#CAG_50 -2.087041 -0.4691949 0.7391099  0.08568188   126377
#CAG_106726 -0.5152433 -0.11912917 0.6958974  0.11787661    18533
#CAG_26 -0.9252039 -0.12923607 0.7220837  0.08233675    16859
#CAG_62 -0.5742917  0.03913569 0.5777953  0.12984145    27507

# The ones obviously linked with high nutrients
agcscore[which(agcscore$RDA2>(0.5)),]
#CAG_112  0.41694778 0.5598896 0.2625683  0.06373306    19164
#CAG_124  0.01666324 0.6506643 0.3284874  0.09524209    26190
#CAG_134 -0.18856372 1.2945381 0.3698938  0.11265678    92550
#CAG_14 -0.27933749 0.6820579 0.4011857  0.08937335    45187
#CAG_176  0.19458562 0.5590770 0.2497336  0.04963427     6733
#CAG_38 -0.43750948 1.0943698 0.3752971  0.08562046    52211

# The ones obviously linked with high temperature
agcscore[which(agcscore$RDA1>(0.5)),]
#CAG_118 0.5027865 -0.4152775 0.3530116   0.1226374     3348
#CAG_12 0.7766346 -0.2250867 0.4498235   0.0999450    59647
#CAG_22 0.6207216 -0.3388147 0.3463599   0.1029897     9528
#CAG_37 0.6317501 -0.6017416 0.3347512   0.1104486     8377
#CAG_42 0.4804619 -0.4427806 0.6463046  0.12400119     9024



# Investigation of other dimensions
   # Axes 1 and 3
scores_13=scoresRDA(res.rda,ax1=1,ax2=3)
samples_13 = scores_13$Sites
samples_13 = as.data.frame(samples_13)
samples_13 = merge(samples_13, Envi_Mat_sel, by = "row.names")
rownames(samples_13) = samples_13[,1]
samples_13 = samples_13[,-1]
samples_13 = merge(samples_13,Group_Var,by="row.names")
rownames(samples_13) = samples_13[,1]
samples_13 = samples_13[,-1]

agcscore_13 = scores_13$Species

enviscore_13=scores_13$Biplot

ggplot() +   geom_hline(yintercept = 0, linetype='dotted') +
  geom_vline(xintercept = 0, linetype='dotted') +
  labs(x = paste0("RDA1 (",
                      round(100 * scores_13$Eigenval.rel[1], 1), "%)"),
           y = paste0("RDA3 (",
                      round(100 * scores_13$Eigenval.rel[3], 1), "%)")) +
  theme(plot.title=element_text(hjust=0.5)) +
  geom_point(aes(x = agcscore_13$RDA1, y = agcscore_13$RDA3),size = 0.5, alpha = 0.4, col = "grey") +
  geom_point(aes(samples_13$RDA1, samples_13$RDA3, col = samples_13$Water_mass_simplified.x,shape=samples_13$MertzGlacier.x), size = 3, alpha = 0.9) +
  geom_segment(data=as.data.frame(enviscore_13),aes(xend = RDA1, yend = RDA3),x=0,y=0,size = 0.5, linetype="F1",color = 'steelblue4',arrow = arrow(length = unit(0.2,"cm"))) +
  geom_text(aes(enviscore_13[,1]*1.05, enviscore_13[,2]*1.05, label=rownames(enviscore_13)), color="steelblue4") +
  scale_colour_manual(values = levels(samples_13$Water_mass_simplified.col)) +
  labs(col = "Water mass") +
  theme_bw()

agcscore_13=merge(agcscore_13,RF,by.x="row.names",by.y="CAG_ID")

ggplot() +   geom_hline(yintercept = 0, linetype='dotted') +
  geom_vline(xintercept = 0, linetype='dotted') +
  labs(x = paste0("RDA1 (",
                      round(100 * scores$Eigenval.rel[1], 1), "%)"),
           y = paste0("RDA3 (",
                      round(100 * scores$Eigenval.rel[3], 1), "%)")) +
  theme(plot.title=element_text(hjust=0.5)) +
  geom_point(aes(x = agcscore_13$RDA1, y = agcscore_13$RDA3, col=agcscore_13$R2_caret,size=agcscore_13$CAG_Size), alpha = 0.6) +
  scale_color_viridis(option="plasma") +
  theme_bw()

# CAGs of interest
agcscore_13[which(agcscore_13$RDA3<(-0.4)),]
#CAG_125 0.4361302 -0.4389633 0.3683815  0.09030785    11038
#CAG_321 0.4650549 -0.4608528 0.2815355  0.09329778     2900
#Alreadya maxis one CAG_37 0.6533071 -0.6321147 0.3347512  0.11044860     8377

agcscore_13[which(agcscore_13$RDA3>(0.4)),]
#CAG_2051 -0.23290963 0.4483729 0.5293704  0.13460756    45205
#CAG_24  0.34127797 0.4208026 0.2768529  0.07547482    15445
#CAG_273  0.38762403 0.5234238 0.5029175  0.08599582     5791
#CAG_44  0.46155708 0.5812630 0.3334760  0.09693514    33274
#CAG_49  0.08126086 0.4080826 0.2155862  0.04263293     7077

   # Axes 1 and 4
scores_14=scoresRDA(res.rda,ax1=1,ax2=4)
samples_14 = scores_14$Sites
samples_14 = as.data.frame(samples_14)
samples_14 = merge(samples_14, Envi_Mat_sel, by = "row.names")
rownames(samples_14) = samples_14[,1]
samples_14 = samples_14[,-1]
samples_14 = merge(samples_14,Group_Var,by="row.names")
rownames(samples_14) = samples_14[,1]
samples_14 = samples_14[,-1]

agcscore_14 = scores_14$Species

enviscore_14=scores_14$Biplot

ggplot() +   geom_hline(yintercept = 0, linetype='dotted') +
  geom_vline(xintercept = 0, linetype='dotted') +
  labs(x = paste0("RDA1 (",
                      round(100 * scores_14$Eigenval.rel[1], 1), "%)"),
           y = paste0("RDA4 (",
                      round(100 * scores_14$Eigenval.rel[4], 1), "%)")) +
  theme(plot.title=element_text(hjust=0.5)) +
  geom_point(aes(x = agcscore_14$RDA1, y = agcscore_14$RDA4),size = 0.5, alpha = 0.4, col = "grey") +
  geom_point(aes(samples_14$RDA1, samples_14$RDA4, col = samples_14$Water_mass_simplified.x,shape=samples_14$MertzGlacier.x), size = 3, alpha = 0.9) +
  geom_segment(data=as.data.frame(enviscore_14),aes(xend = RDA1, yend = RDA4),x=0,y=0,size = 0.5, linetype="F1",color = 'steelblue4',arrow = arrow(length = unit(0.2,"cm"))) +
  geom_text(aes(enviscore_14[,1]*1.05, enviscore_14[,2]*1.05, label=rownames(enviscore_14)), color="steelblue4") +
  scale_colour_manual(values = levels(samples_14$Water_mass_simplified.col)) +
  labs(col = "Water mass") +
  theme_bw()

#Adding R2 info
agcscore_14=merge(agcscore_14,RF,by.x="row.names",by.y="CAG_ID")
ggplot() +   geom_hline(yintercept = 0, linetype='dotted') +
  geom_vline(xintercept = 0, linetype='dotted') +
  labs(x = paste0("RDA1 (",
                      round(100 * scores$Eigenval.rel[1], 1), "%)"),
           y = paste0("RDA4 (",
                      round(100 * scores$Eigenval.rel[4], 1), "%)")) +
  theme(plot.title=element_text(hjust=0.5)) +
  geom_point(aes(x = agcscore_14$RDA1, y = agcscore_14$RDA4, col=agcscore_14$R2_caret,size=agcscore_14$CAG_Size), alpha = 0.6) +
  scale_color_viridis(option="plasma") +
  theme_bw()

agcscore_14[which(agcscore_14$RDA4<(-0.5)),]
#CAG_14 -0.3571304 -0.6481494 0.4011857  0.08937335    45187
#CAG_38 -0.5593519 -0.9946952 0.3752971  0.08562046    52211

agcscore_14[which(agcscore_14$RDA4>0.4),]
#CAG_1069 0.5757562 0.4398006 0.3925455  0.09714608    15704
#CAG_125 0.5391898 0.5314316 0.3683815  0.09030785    11038

   # Axes 1 and 5
scores_15=scoresRDA(res.rda,ax1=1,ax2=5)
samples_15 = scores_15$Sites
samples_15 = as.data.frame(samples_15)
samples_15 = merge(samples_15, Envi_Mat_sel, by = "row.names")
rownames(samples_15) = samples_15[,1]
samples_15 = samples_15[,-1]
samples_15 = merge(samples_15,Group_Var,by="row.names")
rownames(samples_15) = samples_15[,1]
samples_15 = samples_15[,-1]

agcscore_15 = scores_15$Species

enviscore_15=scores_15$Biplot

ggplot() +   geom_hline(yintercept = 0, linetype='dotted') +
  geom_vline(xintercept = 0, linetype='dotted') +
  labs(x = paste0("RDA1 (",
                      round(100 * scores_15$Eigenval.rel[1], 1), "%)"),
           y = paste0("RDA5 (",
                      round(100 * scores_15$Eigenval.rel[5], 1), "%)")) +
  theme(plot.title=element_text(hjust=0.5)) +
  geom_point(aes(x = agcscore_15$RDA1, y = agcscore_15$RDA5),size = 0.5, alpha = 0.4, col = "grey") +
  geom_point(aes(samples_15$RDA1, samples_15$RDA5, col = samples_15$Water_mass_simplified.x,shape=samples_15$MertzGlacier.x), size = 3, alpha = 0.9) +
  geom_segment(data=as.data.frame(enviscore_15),aes(xend = RDA1, yend = RDA5),x=0,y=0,size = 0.5, linetype="F1",color = 'steelblue4',arrow = arrow(length = unit(0.2,"cm"))) +
  geom_text(aes(enviscore_15[,1]*1.05, enviscore_15[,2]*1.05, label=rownames(enviscore_15)), color="steelblue4") +
  scale_colour_manual(values = levels(samples_15$Water_mass_simplified.col)) +
  labs(col = "Water mass") +
  theme_bw()

#Adding R2 info
agcscore_15=merge(agcscore_15,RF,by.x="row.names",by.y="CAG_ID")
ggplot() +   geom_hline(yintercept = 0, linetype='dotted') +
  geom_vline(xintercept = 0, linetype='dotted') +
  labs(x = paste0("RDA1 (",
                      round(100 * scores$Eigenval.rel[1], 1), "%)"),
           y = paste0("RDA5 (",
                      round(100 * scores$Eigenval.rel[5], 1), "%)")) +
  theme(plot.title=element_text(hjust=0.5)) +
  geom_point(aes(x = agcscore_15$RDA1, y = agcscore_15$RDA5, col=agcscore_15$R2_caret,size=agcscore_15$CAG_Size), alpha = 0.6) +
  scale_color_viridis(option="plasma") +
  theme_bw()

#CAGs of interest
agcscore_15[which(agcscore_15$RDA5<(-0.5)),]
#CAG_1069 0.6865443 -0.5611008 0.3925455  0.09714608    15704
#CAG_25 0.3134123 -0.5919113 0.3116166  0.08358431    23763
agcscore_15[which(agcscore_15$RDA5>0.5),]
#CAG_17 0.6075370 0.7243934 0.2791513  0.07667039    34134
#CAG_278 0.4877111 0.9647442 0.2975458  0.13887555     4591
#CAG_53 0.4572232 0.5037414 0.3144350  0.10066276    14968

agcscore[which(agcscore$R2_caret>0.6 & agcscore$CAG_Size>100),]

#CAG_100335 -0.37464336 -0.09006480 0.6743736  0.10307700    14199
#CAG_109240 -0.09135151 -0.02020869 0.6949922  0.09460793     1102
#CAG_29660 -0.03040865 -0.01436784 0.6552114  0.11770181      174
#CAG_42152  0.09755587 -0.08495552 0.6197619  0.10089528      429
#CAG_68443  0.14858934 -0.12368744 0.6014022  0.10898682     1026
#CAG_69320 -0.10794626 -0.02748685 0.6059413  0.13849043     1919
#CAG_75072 -0.25241955 -0.06255918 0.6270375  0.10826231     8987
#CAG_83712 -0.07688978 -0.01401707 0.6135084  0.12950215     1171
#CAG_95210 -0.27326185 -0.06643122 0.7004913  0.08282695     7507


