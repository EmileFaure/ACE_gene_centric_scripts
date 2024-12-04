#Packages
library(data.table)
library(vegan)
library(tidyverse)
library(data.table)
library(viridis)
source("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/R_scripts/RDA/scoresRDA.R")

#Importing Data
#>10% R2
Gene_Mat=fread("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/GM_CAGs_T60MAXGQ_transposed_NZVuniquecut20_FL_RSquared10.tsv",data.table=F)
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
#RsquareAdj(RDA) #76.1%

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
#ggsave("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/01-FIGURES/RDA_Triplot_CAGs_FreeLiving_R2supp10.png",width=29.7,height=21,unit="cm")

#Check R-square value for the clusters in the triplot
RF=fread("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/CAGs_T60MAX_FL_RSquared10_AggregateWithRF.txt",data.table=F)
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
#ggsave("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/01-FIGURES/RDA_CAGplot_FL_R2supp10_colbyR2.png",width=21,height=30,unit="cm")
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
#ggsave("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/01-FIGURES/RDA_CAGplot_FL_R2supp10_onlyEnviSites.png",width=21,height=30,unit="cm")

# Select most informative environmental drivers
set.seed(1994)
#rda.step.both = ordistep(rda(Gene_Mat_hel ~1,data=Envi_Mat_sel),scope=formula(RDA),direction="both",steps=100)
#res.rda = rda(formula=formula(rda.step.both), data=Envi_Mat_sel)
#saveRDS(res.rda,"/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/RDA/RDA_Object_PostVariableSelection_CAGs_FreeLiving_R2supp10.rds")

res.rda=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/RDA/RDA_Object_PostVariableSelection_CAGs_FreeLiving_R2supp10.rds")
RsquareAdj(res.rda) #72.9
#Test significance (Warning, long to run):
#set.seed(1994)
# Full model
#anova.rda <- anova(res.rda) # Significant 0.001
# By axis
#anova.rda.ax <- anova(res.rda, by = "axis", permutations = how(nperm = 499))

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
  geom_text(aes(enviscore[,1]*1.05, enviscore[,2]*1.05, label=rownames(enviscore)), color="steelblue4") +
  scale_colour_manual(values = levels(samples$Water_mass_simplified.col)) +
  labs(col = "Water mass") +
  theme_bw()
#ggsave("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/01-FIGURES/RDA_Triplot_CAGs_FreeLiving_R2supp10_PostVarSel.png",width=23,height=23,unit="cm")

#Adding R2 info
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
ggsave("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/01-FIGURES/RDA_CAGplot_FL_R2supp10_colbyR2_PostVarSel.pdf",width=23,height=23,unit="cm")
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
ggsave("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/01-FIGURES/RDA_CAGplot_FL_R2supp10_onlyEnviSites_PostVarSel.pdf",width=23,height=23,unit="cm")

#Identification of key CAG
# The ones obviously linked with mertz
agcscore[which(agcscore$RDA1<(-1.5)),] 
#CAG_29 -1.680635  1.51727335 0.6683602   0.1132882    38808
#CAG_79 -1.500822  1.69949381 0.7034655   0.1014597    18860

# The ones slightly less Mertz but still
agcscore[which(agcscore$RDA2>0.5),]
#CAG_137 -0.6799384 0.8277440 0.5860565   0.1111682    12036
#CAG_85270 -1.0530279 0.7589338 0.5544879   0.1258320    12097

# The ones linked with deep water masses
agcscore[which(agcscore$RDA1>2),]
#CAG_33 2.568407 0.19697127 0.4016032  0.09331836    41631
#CAG_39 2.962889 0.01269765 0.3852231  0.09242002    30159

# The ones linked with High Temperature / Nitrite
agcscore[which(agcscore$RDA2<(-0.75)),]
# CAG_83  0.01525198 -0.9208024 0.6166796  0.12947590    16318 # More likely teperature
# CAG_35 -0.80720087 -1.3373318 0.3648788  0.12162179    12748
# CAG_34 -0.77884652 -0.9335109 0.4208879  0.09378808    67176
# CAG_22 -1.19926463 -0.7339504 0.4084387  0.11360537    13836 # More likely Nitrite

# The one on axis 1 probably high;y linked to fluorescence
agcscore[which(agcscore$RDA1<(-1)),]
# CAG_131 -1.221765 -0.01165992 0.4549377   0.1133370    11931


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
#Adding R2 info
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
agcscore_13[which(agcscore_13$RDA3<(-1)),]
#CAG_40 -0.4977389 -1.231327 0.3550261  0.09096175   119068
agcscore_13[which(agcscore_13$RDA3>0.4),]
#CAG_83  0.01296282 0.7269998 0.6166796   0.1294759    16318
#CAG_138 -0.04016892 0.4182949 0.5726187   0.1318273     7569

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

# A few samples really stand out on axis 4
samples_14[which(samples_14$RDA4<(-2)),] #Two from 1465 15m

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

#CAGs of interest
agcscore_14[which(agcscore_14$RDA4<(-0.5)),]
#CAG_82 -0.20008039 -1.8085043 0.2152941  0.08915490     5932
#CAG_196 -0.21882677 -1.2286256 0.1945351  0.06850581     8988
#CAG_11  0.07687366 -0.8717204 0.2808033  0.11900107     2763
#CAG_307481 -0.27350874 -0.8515506 0.3476645  0.15681345     3512
#CAG_86 -0.36547548 -0.8086891 0.4469209  0.09660971     2901


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

# Seem to distinguish some deep samples

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
#CAG_67400 2.053048 -0.5893374 0.4070445    0.105873     7062
agcscore_15[which(agcscore_15$RDA5>0.4),]
#CAG_133918  2.0298646 0.7547916 0.3547975  0.07714196    21664
#CAG_558  1.2858515 0.4067378 0.3314069  0.09445123    18018
#+CAG34 which was already standing out on forst two axes


# Finally let's check some CAGs with high R2 that still did not shine
ggplot() +   geom_hline(yintercept = 0, linetype='dotted') +
  geom_vline(xintercept = 0, linetype='dotted') +
  labs(x = paste0("RDA1 (",
                      round(100 * scores$Eigenval.rel[1], 1), "%)"),
           y = paste0("RDA2 (",
                      round(100 * scores$Eigenval.rel[2], 1), "%)")) +
  theme(plot.title=element_text(hjust=0.5)) +
  geom_point(aes(x = RDA1, y = RDA2, col=R2_caret,size=CAG_Size),data=agcscore[which(agcscore$R2_caret>0.6),], alpha = 0.6) +
  scale_color_viridis(option="plasma") +
  theme_bw()

agcscore[which(agcscore$R2_caret>0.6),]
# Many of very small size or singletons --> Not shining on RDA
#CAG_136 -5.636370e-02 -4.237318e-02 0.6292989  0.11988008      925
#CAG_177401 -5.244391e-03 -3.754232e-02 0.6228383  0.12193197      138
#CAG_23 -1.196065e-01 -2.237458e-01 0.6074266  0.11254260     1728
#CAG_237951 -6.837971e-03  3.357411e-03 0.6520745  0.06802662        6
#CAG_244310 -7.445537e-02  6.253446e-02 0.6151278  0.13624729      202
#CAG_330241 -1.224491e-02  1.369168e-02 0.6019417  0.10173199       33
#CAG_38630 -1.277820e-01  1.157820e-01 0.6105032  0.13299861      781
#CAG_73614 -4.032405e-02 -3.936192e-02 0.6651848  0.10668896      249
#CAG_78457 -1.560021e-03 -7.818052e-03 0.6472301  0.10232837        5









