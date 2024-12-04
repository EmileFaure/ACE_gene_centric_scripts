library(tidyverse)

Metadata=readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/metadata_GeneMat_KNN.rds")

Metadata[Metadata$ACE_seq_name=="ACE2_1099_48A_150m",33]="DSW"
Metadata[Metadata$ACE_seq_name=="ACE2_1099_48B_150m",33]="DSW"

Metadata <- mutate(Metadata, MertzGlacier = fct_recode(MertzGlacier, "TRUE"="TRUE.Detroit", "TRUE"="TRUE.Surround"))

saveRDS(Metadata,"/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/Metadata_GeneMat_KNN_MertzTF_1099Corr.rds")
