library(tidyverse)
library(data.table)
source("/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/01_scripts/test_for_enrichment.R")

Annotation=fread("/home/datawork-lmee-intranet-nos/ACE/02-SINGLE-ASSEMBLY/output_anvio_REF/07_GENE_MATRICES/Annotation_Table_AGN_CDH_Tax_KEGG_EGG.tsv", data.table=F,fill=T,sep="\t",header=T) %>%
select(-c(GENE_ID,AGC_ID,AGC_Rep,AGC_Size,AGC_Cat,AGC_Singl_Cat)) %>%
distinct()

Annotation_polar=fread("/home/datawork-lmee-intranet-nos/ACE/10-COMPARE-GENE-TARA/Tara_Mapping_metrics/Annotation_Table_Polar_Genes_CDHit_level.tsv", data.table=F,fill=T,sep="\t",header=F)

names(Annotation_polar) <- names(Annotation)

Annotation <- Annotation[-which(Annotation$CDHit_ID %in% Annotation_polar$CDHit_ID),]

Enrichment_EggNOGCat = test_for_enrichment(Annotation_polar, Annotation, Annotation_Type=best_OG_cat,minimum.occ=1000)
Enrichment_EggNOGDesc = test_for_enrichment(Annotation_polar, Annotation, Annotation_Type=best_OG_desc,minimum.occ=1000)

fwrite(Enrichment_EggNOGCat, "/home/datawork-lmee-intranet-nos/ACE/10-COMPARE-GENE-TARA/Tara_Mapping_metrics/EnrichmentPolar_CDHitLevel_EggnogCAT.tsv",quote=FALSE,sep="\t")
fwrite(Enrichment_EggNOGDesc, "/home/datawork-lmee-intranet-nos/ACE/10-COMPARE-GENE-TARA/Tara_Mapping_metrics/EnrichmentPolar_CDHitLevel_EggNogDesc.tsv",quote=FALSE,sep="\t")

