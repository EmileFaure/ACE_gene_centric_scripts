library(vegan)
library(data.table)
library(ggplot2)
library(tidyverse)
library(reshape2)

#100% Detection

sp100=readRDS("/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/scripts/GENE_MATRIX_FIGURES/AccumVector_100Det_Maximum.RDS")

#90% Detection

sp90=readRDS("/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/scripts/GENE_MATRIX_FIGURES/AccumVector_90Det_Maximum.RDS")

#80% Detection

sp80=readRDS("/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/scripts/GENE_MATRIX_FIGURES/AccumVector_80Det_Maximum.RDS")

#70% Detection

sp70=readRDS("/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/scripts/GENE_MATRIX_FIGURES/AccumVector_70Det_Maximum.RDS")

#60% Detection

sp60=readRDS("/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/scripts/GENE_MATRIX_FIGURES/AccumVector_60Det_Maximum.RDS")

#50% Detection

sp50=readRDS("/home/datawork-lmee-intranet-nos/ACE/Megahit_single_assembly/scripts/GENE_MATRIX_FIGURES/AccumVector_50Det_Maximum.RDS")


dataplot <- data.frame(Sites=sp100$sites, Richness100=sp100$richness, Richness90=sp90$richness,Richness80=sp80$richness, Richness70=sp70$richness, Richness60=sp60$richness, Richness50=sp50$richness)
 
dataplot=melt(dataplot,id.vars=c("Sites"),measured.vars=c("Richness100","Richness90","Richness80","Richness70","Richness60","Richness50"))

AccumGraph = ggplot(data=dataplot) +
  geom_point(aes(x=Sites, y=value, col=variable)) +
  geom_line(aes(x=Sites, y=value, group=variable, col=variable)) +
  scale_color_discrete(labels=c("100%","90%","80%","70%","60%","50%")) +
  labs(y="Number of detected AGNOSTOS Clusters (incl. Singletons)",x="Samples",col="Threshold on detection") +
  theme_bw()

ggsave("/home/datawork-lmee-intranet-nos/ACE/04-GENE-CENTRIC/04-FIGURES/AccumGenes_MaxDet_50to100.pdf", plot = AccumGraph,  width=29.7,height=21,units="cm")

