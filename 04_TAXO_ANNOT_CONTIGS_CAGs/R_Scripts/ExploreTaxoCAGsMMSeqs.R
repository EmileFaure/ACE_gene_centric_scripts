library(tidyverse)

# Import taxonomy, computed with mmseqs taxo on Uniref90 database

Taxo=read.table('/home/datawork-lmee-intranet-nos/ACE/08-TAXO-CAGs-OF-INTEREST/Results_Taxo_Uniref90_138FL_Rfriendly.tsv', sep="\t", header=F,fill=T, quote = '', stringsAsFactors = T)
names(Taxo)=c("Contig","Species","Genus","Family","Order","Class","Phylum","Kingdom","Superkingdom")

# 5 top Class:
TopCla=data.frame(sort(table(Taxo$Class[grepl("^(?!uc_).*$", Taxo$Class, perl = TRUE) & Taxo$Class!="unknown" & Taxo$Class!=""]),decreasing=T)[1:5])
#TopCla=data.frame(sort(table(Taxo$Class),decreasing=T)[1:5])
TopCla$Var1=as.character(TopCla$Var1)
# Add Frequence of the unknown :
TopCla[6,]=c("Unknown",length(Taxo$Class[grepl("^(uc_).*$", Taxo$Class, perl = TRUE)])+length(which(Taxo$Class%in%c("unknown",""))))
TopCla$Freq=as.numeric(TopCla$Freq)
# Add Frequence of the rest
TopCla[7,]=c("Else",nrow(Taxo)-sum(TopCla$Freq))
TopCla$Freq=as.numeric(TopCla$Freq)
TopCla$Var1 = as.factor(TopCla$Var1)

TopCla %>% mutate(Var1=fct_reorder(Var1,desc(Freq))) %>%
  ggplot() + geom_col(aes(x=1, y=Freq, col=Var1, fill=Var1)) +
  labs(fill="Family",col="Family") +
  theme_minimal() +
  theme(text = element_text(size=20))

ggsave("/home/datawork-lmee-intranet-nos/ACE/08-TAXO-CAGs-OF-INTEREST/02_Figures/TaxoBarCla_CAG_138FL.pdf",width=18,height=30,unit="cm")


# 5 top families:
TopFam=data.frame(sort(table(Taxo$Family[grepl("^(?!uc_).*$", Taxo$Family, perl = TRUE) & Taxo$Family!="unknown" & Taxo$Family!=""]),decreasing=T)[1:5])
#TopFam=data.frame(sort(table(Taxo$Family),decreasing=T)[1:5])
TopFam$Var1=as.character(TopFam$Var1)
# Add Frequence of the unknown :
TopFam[6,]=c("Unknown",length(Taxo$Family[grepl("^(uc_).*$", Taxo$Family, perl = TRUE)])+length(which(Taxo$Family%in%c("unknown",""))))
TopFam$Freq=as.numeric(TopFam$Freq)
# Add Frequence of the rest :
TopFam[7,]=c("Else",nrow(Taxo)-sum(TopFam$Freq))
TopFam$Freq=as.numeric(TopFam$Freq)
TopFam$Var1 = as.factor(TopFam$Var1)

TopFam %>% mutate(Var1=fct_reorder(Var1,desc(Freq))) %>%
  ggplot() + geom_col(aes(x=1, y=Freq, col=Var1, fill=Var1)) +
  labs(fill="Family",col="Family") +
  theme_minimal() +
  theme(text = element_text(size=20))

ggsave("/home/datawork-lmee-intranet-nos/ACE/08-TAXO-CAGs-OF-INTEREST/02_Figures/TaxoBarFam_CAG_35FL.pdf",width=18,height=30,unit="cm")

# 5 top genus:
TopGen=data.frame(sort(table(Taxo$Genus[grepl("^(?!uc_).*$", Taxo$Genus, perl = TRUE) & Taxo$Genus!="unknown" & Taxo$Genus!=""]),decreasing=T)[1:5])
#TopGen=data.frame(sort(table(Taxo$Genus),decreasing=T)[1:5])
TopGen$Var1 = as.character(TopGen$Var1)
# Add Frequence of the unknown :
TopGen[6,]=c("Unknown",length(Taxo$Genus[grepl("^(uc_).*$", Taxo$Genus, perl = TRUE)])+length(which(Taxo$Genus%in%c("unknown",""))))
TopGen$Freq=as.numeric(TopGen$Freq)
# Add Frequence of the rest :
TopGen[7,]=c("Else",nrow(Taxo)-sum(TopGen$Freq))
TopGen$Freq=as.numeric(TopGen$Freq)
TopGen$Var1 = as.factor(TopGen$Var1)

TopGen %>% mutate(Var1=fct_reorder(Var1,desc(Freq))) %>%
  ggplot() + geom_col(aes(x=1, y=Freq, col=Var1, fill=Var1)) +
  labs(fill="Genus",col="Genus") +
  theme_minimal() +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  theme(text = element_text(size=20))

ggsave("/home/datawork-lmee-intranet-nos/ACE/08-TAXO-CAGs-OF-INTEREST/02_Figures/TaxoBarGen_CAG_35FL.pdf",width=18,height=30,unit="cm")



