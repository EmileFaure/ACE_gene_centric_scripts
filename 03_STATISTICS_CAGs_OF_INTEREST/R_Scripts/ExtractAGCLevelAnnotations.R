library(tidyverse)
library(data.table)

Annotation=fread("/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/Annotation_Table_AGN_CDH_Tax_KEGG_EGG.tsv", data.table=F,fill=T,sep="\t",header=T)

Mode <- function(x) {
  ux <- na.omit(unique(x))
  ux[which.max(tabulate(match(x, ux)))]
} # Warning, if two modes, only the first appearing will be kept

Annotation = Annotation %>% select(-c("GENE_ID","DOMAIN","CDHit_ID")) %>%
 group_by(AGC_ID) %>%
 summarise_all(Mode)

fwrite(Annotation, "/home/datawork-lmee-intranet-nos/ACE/07-STATS-CAGs-OF-INTEREST/Annotation_Table_AGNLevel",quote=FALSE,sep="\t")
