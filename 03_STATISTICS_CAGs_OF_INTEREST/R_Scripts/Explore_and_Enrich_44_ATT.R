library(tidyverse)
library(data.table)
library(vegan)
library(reshape2)
library(treemap)
source("/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/01_scripts/test_for_enrichment.R")

# Import data
Gene_Mat=fread("/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/GM_AGN_T60MAXGQ_transposed_CAG_44_ATT.tsv", data.table=F)
names(Gene_Mat)[1]="Sample_ID"
Annotation=fread("/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/Annotation_Table_AGN_CDH_Tax_KEGG_EGG_CAG_44ATT.tsv", data.table=F,fill=T,sep="\t",header=T)
Envi_Mat=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/Metadata_GeneMat_KNN_MertzTF_1099Corr.rds")
Group_Var=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/metadata_metaG_grouping_variables.rds")
RF=fread("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/RF_AGC_NZVuniquecut20_T60MAX_ATT_All.txt",data.table=F)
RF=RF[which(RF$V1 %in% names(Gene_Mat)),]

# Exploration of RF importance
rank.rf.imp = as.data.frame(t(apply(as.matrix(RF[,c(6:ncol(RF))]), 1, rank)))
rank.rf.imp = 50-rank.rf.imp
boxplot = reshape2::melt(as.matrix(rank.rf.imp))
boxplot$Var2=sub("Importance.","",boxplot$Var2)
# All variables
#bpall = ggplot() + geom_boxplot(aes(x = reorder(Var2,value,FUN=median), y = value), data = boxplot) +
#  theme(axis.text.x = element_text(angle = 90, hjust = 1, face="bold")) +
#  labs(x = "", y = "Rank of importance in random forest regressions")
#bpall
#ggsave("/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/02-FIGURES/Importance_Boxplot_CAGXXXX.png",width=29.7,height=21,unit="cm")
# Top Ten Only
boxplotTop <- boxplot %>% group_by(Var2) %>%
 summarise(median=median(value)) %>%
 arrange(median) %>%
 pull(Var2)
boxplotTop=boxplotTop[1:10]
bptop=ggplot() + geom_boxplot(aes(x = reorder(Var2,value,FUN=median), y = value), data = boxplot[which(boxplot$Var2 %in% boxplotTop),]) +
  scale_y_reverse() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face="bold"),text = element_text(size=20)) +
  labs(x = "", y = "Rank of importance in random forest regressions")
bptop
ggsave("/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/02-FIGURES/Importance_Boxplot_CAG_44ATT_Top10.pdf",width=18,height=25,unit="cm")

# Exploration of abundance at CAG level in response to environment
var_envi1=sym(boxplotTop[1])
var_envi2=sym(boxplotTop[2])
CAG_Abund=data.frame(cbind(Gene_Mat[,1],rowMeans(Gene_Mat[,-1])))
names(CAG_Abund)=c("Sample_ID","Abundance")
CAG_Abund$Sample_ID=as.factor(CAG_Abund$Sample_ID)
CAG_Abund$Abundance=as.numeric(CAG_Abund$Abundance)
CAG_Abund=merge(CAG_Abund,Envi_Mat,by.x="Sample_ID", by.y="ACE_seq_name")
ggplot() + geom_point(aes(y=Abundance, x=!!var_envi1, shape=MertzGlacier, col=!!var_envi2), size=4,alpha=0.8, data=CAG_Abund) +
 theme_bw() +
 theme(text = element_text(size=20))
ggsave("/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/02-FIGURES/AbundEnviPlot_CAG_44ATT.png",width=23,height=21,unit="cm")

####  Exploration of Annotations at AGC level
AGC_level_annot=select(Annotation,c("AGC_ID","AGC_Size","AGC_Cat","AGC_Singl_Cat"))
AGC_level_annot=AGC_level_annot[!duplicated(AGC_level_annot),]
Known=as.numeric(sum(table(AGC_level_annot$AGC_Cat)["K"],table(AGC_level_annot$AGC_Singl_Cat)["K"],na.rm=T))
Known_With_PFAM=as.numeric(sum(table(AGC_level_annot$AGC_Cat)["KWP"],table(AGC_level_annot$AGC_Singl_Cat)["KWP"],na.rm=T))
Genomic_Unknown=as.numeric(sum(table(AGC_level_annot$AGC_Cat)["GU"],table(AGC_level_annot$AGC_Singl_Cat)["GU"],na.rm=T))
Environmental_Unknown=as.numeric(sum(table(AGC_level_annot$AGC_Cat)["EU"],table(AGC_level_annot$AGC_Singl_Cat)["EU"],na.rm=T))
Total=Known+Known_With_PFAM+Genomic_Unknown+Environmental_Unknown
Tot_Known=Known+Known_With_PFAM
Tot_Unknown=Genomic_Unknown+Environmental_Unknown
# Treemap representation
agnostos=data.frame(Type=c(paste0("Known (",round(Tot_Known/Total*100,1),"%)"),paste0("Known (",round(Tot_Known/Total*100,1),"%)"),paste0("Unknown (",round(Tot_Unknown/Total*100,1),"%)"),paste0("Unknown (",round(Tot_Unknown/Total*100,1),"%)")),
Cluster_Type=c(paste0("Known\n(",round(Known/Total*100,1),"%)"),paste0("Known with PFAM\n(",round(Known_With_PFAM/Total*100,1),"%)"),paste0("Genomic Unknown\n(",round(Genomic_Unknown/Total*100,1),"%)"),paste0("Environmental Unknown\n(",round(Environmental_Unknown/Total*100,1),"%)")),
                    Clusters=c(Known,Known_With_PFAM,Genomic_Unknown,Environmental_Unknown),colors=c("#339966","#8fbc8f","#3399CC","#FF6666"))

png(filename="/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/02-FIGURES/CAG_44ATT_AGNOSTOSTreeMap.png",width=20, height=20,units="cm",res=1200)
treemap(agnostos, index=c("Type","Cluster_Type"),     vSize="Clusters", type="color",
        fontsize.labels=c(20,15),                # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
        fontcolor.labels=c("black","white"),    # Color of labels
        fontface.labels=c(2,4),                  # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
        bg.labels = 0,
        align.labels=list(
          c("left", "top"),
          c("center", "center")
        ),                                   # Where to place labels in the rectangle?
        overlap.labels=0.5,                      # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
        inflate.labels=F,                        # If true, labels are bigger when rectangle is bigger.
        border.col=c("black","white"),             # Color of borders of groups, of subgroups, of subsubgroups ....
        border.lwds=c(7,2),                 # Width of colors
        vColor="colors",
        title = ""
)
dev.off()

## Cheking MAGs taxonomy
  # In CAG 29 there are 3,283,742 contigs. 552,722 are found in MAGs, among which 156,206 are in flavobacterriaceae MAGs
MAG_Tax=fread("/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/Taxo_MAGs_CAG_44ATT.tsv", data.table=F, header=F, fill=T, stringsAsFactors=T)
names(MAG_Tax)=c("ContigID","MAG_ID","REP_MAG_ID","Domain","Phylum","Class","Order","Family","Genus","Species")
list_contigs=fread("/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/list_contigs_CAG_44ATT.tsv",data.table=F, header=F, fill=T, stringsAsFactors=T)
ncontigs=nrow(list_contigs)
rm(list_contigs)

MAG_plot_Tot=data.frame(
Classif=as.factor(c("In_MAG","Not_In_MAG")),
Number=c(nrow(MAG_Tax),ncontigs-nrow(MAG_Tax)))

ggplot() + geom_col(aes(x=1, y=Number, col=Classif, fill=Classif), data=MAG_plot_Tot) +
 scale_fill_manual(values=c("steelblue","indianred")) +
 scale_color_manual(values=c("steelblue","indianred")) +
 theme_minimal() +
 theme(text = element_text(size=20))

ggsave("/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/02-FIGURES/TaxoBarTot_CAG_44ATT.pdf",width=18,height=30,unit="cm")

# 5 top families:
TopFam=data.frame(sort(table(MAG_Tax$Family),decreasing=T)[1:5])
# Add Frequence of the rest :
TopFam$Var1=sub("f__","",TopFam$Var1)
TopFam[6,]=c("Else",nrow(MAG_Tax)-sum(data.frame(sort(table(MAG_Tax$Family),decreasing=T)[1:5])$Freq))
TopFam$Freq=as.numeric(TopFam$Freq)
TopFam$Var1 = as.factor(TopFam$Var1)

TopFam %>% mutate(Var1=fct_reorder(Var1,desc(Freq))) %>%
ggplot() + geom_col(aes(x=1, y=Freq, col=Var1, fill=Var1)) +
 labs(fill="Family",col="Family") +
 theme_minimal() +
 theme(text = element_text(size=20))

ggsave("/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/02-FIGURES/TaxoBarFam_CAG_44ATT.pdf",width=18,height=30,unit="cm")

# 5 top genus:
TopGen=data.frame(sort(table(MAG_Tax$Genus),decreasing=T)[1:5])
# Add Frequence of the rest :
TopGen$Var1=sub("g__","",TopGen$Var1)
TopGen[6,]=c("Else",nrow(MAG_Tax)-sum(data.frame(sort(table(MAG_Tax$Genus),decreasing=T)[1:5])$Freq))
TopGen$Freq=as.numeric(TopGen$Freq)
TopGen$Var1 = as.factor(TopGen$Var1)

TopGen %>% mutate(Var1=fct_reorder(Var1,desc(Freq))) %>%
ggplot() + geom_col(aes(x=1, y=Freq, col=Var1, fill=Var1)) +
 labs(fill="Genus",col="Genus") +
 theme_minimal() +
 scale_color_brewer(palette="Set1") +
 scale_fill_brewer(palette="Set1") +
 theme(text = element_text(size=20))
ggsave("/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/02-FIGURES/TaxoBarGen_CAG_44ATT.pdf",width=18,height=30,unit="cm")


###  FUNCTIONAL ENRICHMENT OF THE CAG  ###
# Here done at AGC level, can also be achieved at CD-Hit level

Annotation=fread("/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/Annotation_Table_AGNLevel.tsv", data.table=F,fill=T,sep="\t",header=T)
CAG_AGCList=fread("/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/CAG_44_ATTRsquared15.txt", data.table=F,fill=T,sep="\t",header=F)
names(CAG_AGCList)=c("CAG_ID","AGC_ID")

Annotation_CAG=Annotation[which(Annotation$AGC_ID %in% CAG_AGCList$AGC_ID),]
Annotation=Annotation[-which(Annotation$AGC_ID %in% CAG_AGCList$AGC_ID),]

# Check functional enrichment
# a- EggNOG Cat

Enrichment_EggNOGCat = test_for_enrichment(Annotation_CAG, Annotation, Annotation_Type=best_OG_cat,minimum.occ=1000)

fwrite(Enrichment_EggNOGCat, "/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/EnrichmentCAG_44ATT_AGCLevel.tsv",quote=FALSE,sep="\t")

Enrichment_EggNOGCat[which(Enrichment_EggNOGCat$Enriched==T),] %>%
mutate(Annotation = fct_reorder(Annotation, desc(Corrected_p_value))) %>%
filter(Annotation != "-" & Annotation != "S" & Annotation != "") %>% #Removes unknown (should be dealt with independantly as multiple annotations indicate unknown
arrange(Corrected_p_value) %>%
slice(1:10) %>%
ggplot() + geom_col(aes(y=Annotation,x=Odds_Ratio)) +
 theme_bw() +
 theme(text = element_text(size=20))
ggsave("/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/02-FIGURES/EnrichmentCOGCat_MinOcc1000_AGCLevel_CAG_44ATT_Top10Only.pdf",width=18,height=30,unit="cm")

# b- EggNOG Desc

Enrichment_EggNOGDesc = test_for_enrichment(Annotation_CAG, Annotation, Annotation_Type=best_OG_desc,minimum.occ=1000)

fwrite(Enrichment_EggNOGDesc, "/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/EnrichmentCAG_44ATT_AGCLevel_EggNogDesc.tsv",quote=FALSE,sep="\t")

if (nrow(Enrichment_EggNOGDesc[which(Enrichment_EggNOGDesc$Enriched==T),])!=0) {
Enrichment_EggNOGDesc[which(Enrichment_EggNOGDesc$Enriched==T),] %>%
mutate(Annotation = fct_reorder(Annotation, desc(Corrected_p_value))) %>%
filter(Annotation != "-" & Annotation != "") %>%
arrange(Corrected_p_value) %>%
slice(1:10) %>%
ggplot() + geom_col(aes(y=Annotation,x=Odds_Ratio)) +
 theme_bw() +
 scale_y_discrete(labels=function(x) str_wrap(x,width=16)) +
 theme(text = element_text(size=12))
ggsave("/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/02-FIGURES/EnrichmentCOGDesc_MinOcc1000_AGCLevel_CAG_44ATT_Top10Only.pdf",width=18,height=30,unit="cm")
}

# c- KEGG KO

Enrichment_KEGG = test_for_enrichment(Annotation_CAG, Annotation, Annotation_Type=KEGG_KO,minimum.occ=100)

fwrite(Enrichment_KEGG, "/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/EnrichmentCAG_44ATT_AGCLevel_KeggKO.tsv",quote=FALSE,sep="\t")

Enrichment_KEGG[which(Enrichment_KEGG$Enriched==T),] %>%
mutate(Annotation = fct_reorder(Annotation, desc(Corrected_p_value))) %>%
arrange(Corrected_p_value) %>%
slice(1:10) %>%
ggplot() + geom_col(aes(y=Annotation,x=Odds_Ratio)) +
 theme_bw() +
 theme(text = element_text(size=15))
ggsave("/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/02-FIGURES/EnrichmentKEGG_MinOcc100_AGCLevel_CAG_44ATT_Top10Only.pdf",width=18,height=30,unit="cm")

