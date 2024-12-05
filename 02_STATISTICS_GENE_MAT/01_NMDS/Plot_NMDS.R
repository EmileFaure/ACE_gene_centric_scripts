library(ggplot2)
library(vegan)
library(tidyverse)
library(viridis)

MetaData = readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/Metadata_GeneMat_KNN_MertzTF_1099Corr.rds")
NMDSData = readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/NMDS_AGC_NZVuniquecut20.rds")
ColorData = readRDS("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/metadata_metaG_grouping_variables.rds")

df.NMDS <- scores(NMDSData, display = "sites") %>% as.data.frame()

df.NMDS <- merge(df.NMDS,MetaData,by.x="row.names",by.y="ACE_seq_name")

df.NMDS$Size_fraction <- factor(df.NMDS$Size_fraction, levels = c("0.2-3 µm", "0.2-40 µm", ">3 µm"))

plot.NMDS.wm <- ggplot(data = df.NMDS,
                    aes(x = NMDS1, y = NMDS2, col = Water_mass_simplified, shape = MertzGlacier, size = Size_fraction)) +
  geom_point() +
  scale_colour_manual(values = levels(ColorData$Water_mass_simplified.col)) +
  labs(col = "Water Mass", shape = "Mertz Glacier Sample") +
  theme_bw() +
  theme(legend.background = element_blank(),
        legend.box.background= element_rect(colour="black"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) +
  guides(shape = guide_legend(ncol=2))
plot.NMDS.wm
ggsave("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/Figures/plot_NMDS_wm.pdf")

plot.NMDS.lh <- ggplot(data = df.NMDS,
                    aes(x = NMDS1, y = NMDS2, col = Longhurst_Prov, shape = MertzGlacier, size = Size_fraction)) +
  geom_point() +
  labs(col = "Lonhurst biogeographical provinces", shape = "Mertz Glacier Sample") +
  theme_bw() +
  theme(legend.background = element_blank(),
        legend.box.background= element_rect(colour="black"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) +
  guides(shape = guide_legend(ncol=2))
plot.NMDS.lh
ggsave("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/Figures/plot_NMDS_lh.pdf")

plot.NMDS.temp <- ggplot(data = df.NMDS,
                    aes(x = NMDS1, y = NMDS2, col = Temperature_ITS90, shape = MertzGlacier, size = Size_fraction)) +
  geom_point() +
  labs(col = "Temperature", shape = "Mertz Glacier Sample") +
  scale_colour_viridis() +
  theme_bw() +
  theme(legend.background = element_blank(),
        legend.box.background= element_rect(colour="black"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) +
  guides(shape = guide_legend(ncol=2))
plot.NMDS.temp
ggsave("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/Figures/plot_NMDS_temp.pdf")


plot.NMDS.plot.NMDS.lv8 <- ggplot(data = df.NMDS,
                    aes(x = NMDS1, y = NMDS2, col = LV8_mean, shape = MertzGlacier, size = Size_fraction)) +
  geom_point() +
  labs(col = "LV8_mean", shape = "Mertz Glacier Sample") +
  scale_colour_viridis() +
  theme_bw() +
  theme(legend.background = element_blank(),
        legend.box.background= element_rect(colour="black"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) +
  guides(shape = guide_legend(ncol=2))
plot.NMDS.plot.NMDS.lv8

plot.NMDS.plot.NMDS.lv6 <- ggplot(data = df.NMDS,
                    aes(x = NMDS1, y = NMDS2, col = LV6_mean, shape = MertzGlacier, size = Size_fraction)) +
  geom_point() +
  labs(col = "LV6_mean", shape = "Mertz Glacier Sample") +
  scale_colour_viridis() +
  theme_bw() +
  theme(legend.background = element_blank(),
        legend.box.background= element_rect(colour="black"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) +
  guides(shape = guide_legend(ncol=2))
plot.NMDS.plot.NMDS.lv6

plot.NMDS.plot.NMDS.Fe <- ggplot(data = df.NMDS,
                    aes(x = NMDS1, y = NMDS2, col = dFe_nmol.kg, shape = MertzGlacier, size = Size_fraction)) +
  geom_point() +
  labs(col = "dFe_nmol.kg", shape = "Mertz Glacier Sample") +
  scale_colour_viridis() +
  theme_bw() +
  theme(legend.background = element_blank(),
        legend.box.background= element_rect(colour="black"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) +
  guides(shape = guide_legend(ncol=2))
plot.NMDS.plot.NMDS.Fe

plot.NMDS.plot.NMDS.d13C <- ggplot(data = df.NMDS,
                    aes(x = NMDS1, y = NMDS2, col = d13C, shape = MertzGlacier, size = Size_fraction)) +
  geom_point() +
  labs(col = "LV1_mean", shape = "Mertz Glacier Sample") +
  scale_colour_viridis() +
  theme_bw() +
  theme(legend.background = element_blank(),
        legend.box.background= element_rect(colour="black"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) +
  guides(shape = guide_legend(ncol=2))
plot.NMDS.plot.NMDS.d13C
ggsave("/home/datawork-lmee-intranet-nos/ACE/06-STATS-GENE-MAT/Figures/plot_NMDS_d13C.pdf")




