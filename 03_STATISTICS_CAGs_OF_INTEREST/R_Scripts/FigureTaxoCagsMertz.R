library(tidyverse)

taxo_cag79 <- data.frame(Contigs_Number=c(349870,
                                          302760,
                                          15777,
                                          17185,
                                          25269,
                                          249772),
                         Family=c(
                           "Flavobacteriaceae",
                           "Unknown",
                           "Cryomorphaceae",
                           "Roseobacteraceae",
                           "Crocinitomicaceae",
                           "Else"
                         ),
                         CAG=rep("CAG 79",6))

taxo_cag29 <- data.frame(Contigs_Number=c(219592,
                                          209269,
                                          172653,
                                          101260,
                                          858771,
                                          502537),
                         Family=c(
                           "Porticoccaceae",
                           "Oceanospirillaceae",
                           "Flavobacteriaceae",
                           "Rhodobacteraceae",
                           "Unknown",
                           "Else"
                         ),
                         CAG=rep("CAG 29",6))

taxo_cag137 <- data.frame(Contigs_Number=c(86686,
                                          52483,
                                          63301,
                                          41577,
                                          3393,
                                          83695),
                         Family=c(
                           "Roseobacteraceae",
                           "Rhodobacteraceae",
                           "Unknown",
                           "Flavobacteriaceae",
                           "Crocinitomicaceae",
                           "Else"
                         ),
                         CAG=rep("CAG 137",6))

taxo_cag85270 <- data.frame(Contigs_Number=c(520183,
                                             264918,
                                             122077,
                                             75543,
                                             40689,
                                           492590),
                          Family=c(
                            "Unknown",
                            "Flavobacteriaceae",
                            "Rhodobacteraceae",
                            "Roseobacteraceae",
                            "Crocinitomicaceae",
                            "Else"
                          ),
                          CAG=rep("CAG 85270",6))

taxo_fam <- bind_rows(taxo_cag79,taxo_cag29,taxo_cag137,taxo_cag85270)

taxo_fam %>% group_by(CAG) %>%
  mutate(ContigsPerc=Contigs_Number/sum(Contigs_Number)) %>%
  ungroup() %>%
  mutate(Family=factor(Family, levels=c("Unknown","Else","Flavobacteriaceae","Rhodobacteraceae",
                                        "Roseobacteraceae", "Crocinitomicaceae", "Porticoccaceae",
                                        "Oceanospirillaceae","Cryomorphaceae"))) %>%
  mutate(CAG=factor(CAG,levels=c("CAG 29","CAG 79","CAG 137", "CAG 85270"))) %>%
  ggplot() + geom_col(aes(x=CAG, y=ContigsPerc, col=Family, fill=Family)) +
  labs(fill="Family",col="Family") +
  scale_fill_manual(values=c("black","grey","indianred","darkorchid4","darkorchid2","lightblue","forestgreen","blue","orange")) +
  scale_color_manual(values=c("black","grey","indianred","darkorchid4","darkorchid2","lightblue","forestgreen","blue","orange")) +
  theme_minimal() +
  theme(text = element_text(size=20))

### KRAKEN GTDB

taxo_cag79 <- data.frame(Contigs_Number=c(578891,
                                          43902,
                                          39238,
                                          35348,
                                          19150,
                                          244104),
                         Family=c(
                           "Flavobacteriaceae",
                           "Unknown",
                           "BACL11",
                           "Rhodobacteraceae",
                           "UA16",
                           "Else"
                         ),
                         CAG=rep("CAG 79",6))

taxo_cag29 <- data.frame(Contigs_Number=c(367282,
                                          275654,
                                          265275,
                                          189611,
                                          143704,
                                          1072328),
                         Family=c(
                           "Flavobacteriaceae",
                           "Porticoccaceae",
                           "Nitrincolaceae",
                           "Rhodobacteraceae",
                           "Unknown",
                           "Else"
                         ),
                         CAG=rep("CAG 29",6))

taxo_fam <- bind_rows(taxo_cag79,taxo_cag29)

taxo_fam %>% group_by(CAG) %>%
  mutate(ContigsPerc=Contigs_Number/sum(Contigs_Number)) %>%
  ungroup() %>%
  mutate(Family=factor(Family, levels=c("Unknown","Else","Flavobacteriaceae","Rhodobacteraceae",
                                        "Nitrincolaceae", "BACL11", "Porticoccaceae",
                                        "UA16"))) %>%
  mutate(CAG=factor(CAG,levels=c("CAG 29","CAG 79"))) %>%
  ggplot() + geom_col(aes(x=CAG, y=ContigsPerc, col=Family, fill=Family)) +
  labs(fill="Family",col="Family") +
  scale_fill_manual(values=c("black","grey","indianred","darkorchid4","lightblue","blue","forestgreen","orange")) +
  scale_color_manual(values=c("black","grey","indianred","darkorchid4","lightblue","blue","forestgreen","orange")) +
  theme_minimal() +
  theme(text = element_text(size=20))
