library(vegan)
library(data.table)
library(ggplot2)
library(tidyverse)

#100% Detection

#GeneMat=fread("/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES/AGNOSTOS_GQ_CLSTRLVL_GENE_MAT_COV_MAXIMUMDETECTHRESH_100percent.tsv", sep="\t", header=TRUE,stringsAsFactors=TRUE)
#GeneMat=as.data.frame(GeneMat)
#row.names(GeneMat) = GeneMat[,1]
#GeneMat=GeneMat[,-1]
#GeneMat=as.data.frame(t(GeneMat))
# Collector vector computation
#sp100 <- specaccum(GeneMat,"collector")
#rm(GeneMat)
#saveRDS(sp100,file="/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/scripts/GENE_MATRIX_FIGURES/AccumVector_100Det_Maximum.RDS")

#90% Detection

#GeneMat=fread("/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES/AGNOSTOS_GQ_CLSTRLVL_GENE_MAT_COV_MAXIMUMDETECTHRESH_90percent.tsv", sep="\t", header=TRUE,stringsAsFactors=TRUE)
#GeneMat=as.data.frame(GeneMat)
#row.names(GeneMat) = GeneMat[,1]
#GeneMat=GeneMat[,-1]
#GeneMat=as.data.frame(t(GeneMat))
# Collector vector computation
#sp90 <- specaccum(GeneMat,"collector")
#rm(GeneMat)
#saveRDS(sp90,file="/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/scripts/GENE_MATRIX_FIGURES/AccumVector_90Det_Maximum.RDS")

#80% Detection

GeneMat=fread("/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES/AGNOSTOS_GQ_CLSTRLVL_GENE_MAT_COV_MAXIMUMDETECTHRESH_80percent.tsv", sep="\t", header=TRUE,stringsAsFactors=TRUE)
GeneMat=as.data.frame(GeneMat)
row.names(GeneMat) = GeneMat[,1]
GeneMat=GeneMat[,-1]
GeneMat=as.data.frame(t(GeneMat))
# Collector vector computation
sp80 <- specaccum(GeneMat,"collector")
rm(GeneMat)
saveRDS(sp80,file="/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/scripts/GENE_MATRIX_FIGURES/AccumVector_80Det_Maximum.RDS")

#70% Detection

GeneMat=fread("/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES/AGNOSTOS_GQ_CLSTRLVL_GENE_MAT_COV_MAXIMUMDETECTHRESH_70percent.tsv", sep="\t", header=TRUE,stringsAsFactors=TRUE)
GeneMat=as.data.frame(GeneMat)
row.names(GeneMat) = GeneMat[,1]
GeneMat=GeneMat[,-1]
GeneMat=as.data.frame(t(GeneMat))
# Collector vector computation
sp70 <- specaccum(GeneMat,"collector")
rm(GeneMat)
saveRDS(sp70,file="/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/scripts/GENE_MATRIX_FIGURES/AccumVector_70Det_Maximum.RDS")

#60% Detection

GeneMat=fread("/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES/AGNOSTOS_GQ_CLSTRLVL_GENE_MAT_COV_MAXIMUMDETECTHRESH_60percent.tsv", sep="\t", header=TRUE,stringsAsFactors=TRUE)
GeneMat=as.data.frame(GeneMat)
row.names(GeneMat) = GeneMat[,1]
GeneMat=GeneMat[,-1]
GeneMat=as.data.frame(t(GeneMat))
# Collector vector computation
sp60 <- specaccum(GeneMat,"collector")
rm(GeneMat)
saveRDS(sp60,file="/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/scripts/GENE_MATRIX_FIGURES/AccumVector_60Det_Maximum.RDS")

#50% Detection

GeneMat=fread("/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/output_anvio_REF/07_GENE_MATRICES/AGNOSTOS_GQ_CLSTRLVL_GENE_MAT_COV_MAXIMUMDETECTHRESH_50percent.tsv", sep="\t", header=TRUE,stringsAsFactors=TRUE)
GeneMat=as.data.frame(GeneMat)
row.names(GeneMat) = GeneMat[,1]
GeneMat=GeneMat[,-1]
GeneMat=as.data.frame(t(GeneMat))
# Collector vector computation
sp50 <- specaccum(GeneMat,"collector")
rm(GeneMat)
saveRDS(sp50,file="/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/scripts/GENE_MATRIX_FIGURES/AccumVector_50Det_Maximum.RDS")


dataplot <- data.frame(Sites=sp100$sites, Richness100=sp100$richness, Richness90=sp90$richness,Richness80=sp80$richness, Richness70=sp70$richness, Richness60=sp60$richness, Richness50=sp50$richness)
 
dataplot=melt(dataplot,id.vars=c("Sites"),measured.vars=c("Richness100","Richness90","Richness80","Richness70","Richness60","Richness50"))

AccumGraph = ggplot(data=dataplot) +
  geom_point(aes(x=Sites, y=value, col=variable)) +
  geom_line(aes(x=Sites, y=Richness, group=variable, col=variable)) +
  labs(y="Number of detected AGNOSTOS Clusters (incl. Singletons)",x="Samples",col="Threshold on detection") +
  theme_bw()

ggsave("/home/datawork-lmee-intranet-nos/ACE/04-GENE-CENTRIC/04-FIGURES/AccumGenes_MaxDet_50to100.pdf", plot = AccumGraph,  width=29.7,height=21,units="cm")

