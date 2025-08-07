require(tidyverse)
require(vegan)
library(corrplot)
library(ComplexHeatmap)
library(circlize)  # for color functions
library(viridis)
library(RColorBrewer)
library(data.table)
source("/Users/emifaure/Documents/TONGA/StageMilo/PetB_Corr_Heatmaps_RDA/scoresRDA.R")

########## Load metadata ############

meta <- readRDS("/Users/emifaure/Documents/ACE/Metadata_GeneMat_KNN_MertzTF_1099Corr.rds")
Group_Var=readRDS("/Users/emifaure/Documents/ACE/metadata_metaG_grouping_variables.rds")

########## Load DMSP abundance ############

DMSP <- read.table("/Users/emifaure/Documents/ACE/DMSP/GM_AGN_TMAX60GQ_transposed_NZVuniquecut20_FreeLiving_DMSPLyase.tsv", header = TRUE, sep="\t")

# Select metadata from samples :
meta_DMSP <- meta[which(meta$ACE_seq_name %in% row.names(DMSP)),]
rownames(meta_DMSP) <- meta_DMSP$ACE_seq_name
meta_DMSP <- meta_DMSP[,-1]
# Match row order
meta_DMSP <- meta_DMSP[match(row.names(DMSP), row.names(meta_DMSP)),]
#extract numeric for correlation computations :
meta_DMSP_quanti <- meta_DMSP %>% select(where(is.numeric))
#The sign of some latent variables needs to be inversed to be interpreted as in the Landwehr et al. study
#Mail from Landwehr : "We have chosen to change the signs of some of the LVs that where produced by the algorithm, in order to make them align with certain physical processes"
#Also note that here LV are from 0 to 13, like in the data provided by Landwehr et al., need to add one to each LV number to match their study.
meta_DMSP_quanti$LV3_mean = -meta_DMSP_quanti$LV3_mean
meta_DMSP_quanti$LV12_mean = -meta_DMSP_quanti$LV12_mean
meta_DMSP_quanti$LV4_mean = -meta_DMSP_quanti$LV4_mean
meta_DMSP_quanti$LV11_mean = -meta_DMSP_quanti$LV11_mean
meta_DMSP_quanti$LV10_mean = -meta_DMSP_quanti$LV10_mean
meta_DMSP_quanti$LV5_mean = -meta_DMSP_quanti$LV5_mean
meta_DMSP_quanti$LV9_mean = -meta_DMSP_quanti$LV9_mean

########## Explore DMSP abundance ############

# Correlation-based figures
DMSP %>% rownames_to_column("sample") %>% 
  pivot_longer(names_to = "AGC", cols = -sample) %>%
  ggplot() +
  geom_violin(aes(x=sample, y=value))

DMSP_Sum <- as.data.frame(rowSums(DMSP))

# Initialize vectors to store results
cor_values <- numeric(ncol(meta_DMSP_quanti))
p_values <- numeric(ncol(meta_DMSP_quanti))

# Loop over each environmental variable (column)
for (i in seq_along(meta_DMSP_quanti)) {
  test_result <- cor.test(DMSP_Sum[,1], meta_DMSP_quanti[, i], method = "pearson")
  cor_values[i] <- test_result$estimate
  p_values[i] <- test_result$p.value
}

# Adjust p-values using Benjamini-Hochberg (FDR)
p_adj_bh <- p.adjust(p_values, method = "BH")

# Create a summary data frame
results <- data.frame(
  Variable = colnames(meta_DMSP_quanti),
  Correlation = cor_values,
  P_value = p_values,
  P_adj_BH = p_adj_bh
)

results[results$P_adj_BH<0.05,]

# Same but by AGC :
# Initialize matrices
n_abund <- ncol(DMSP)
n_env <- ncol(meta_DMSP_quanti)

cor_matrix <- matrix(NA, nrow = n_abund, ncol = n_env)
p_matrix <- matrix(NA, nrow = n_abund, ncol = n_env)

# Loop over all combinations of abundance vs environment
for (i in 1:n_abund) {
  for (j in 1:n_env) {
    test <- cor.test(DMSP[, i], meta_DMSP_quanti[, j], method = "pearson")
    cor_matrix[i, j] <- test$estimate
    p_matrix[i, j] <- test$p.value
  }
}

# Assign row/column names for clarity
rownames(cor_matrix) <- colnames(DMSP)
colnames(cor_matrix) <- colnames(meta_DMSP_quanti)

rownames(p_matrix) <- colnames(DMSP)
colnames(p_matrix) <- colnames(meta_DMSP_quanti)

# Adjust
p_adj_vector <- p.adjust(as.vector(p_matrix), method = "BH")
p_adj_matrix <- matrix(p_adj_vector, nrow = n_abund, ncol = n_env)
rownames(p_adj_matrix) <- colnames(DMSP)
colnames(p_adj_matrix) <- colnames(meta_DMSP_quanti)

# Define color scale
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Create matrix of significance labels (e.g., "*" for p < 0.05)
sig_labels <- ifelse(p_adj_matrix > 0.05, "x", "")

Heatmap(
  cor_matrix,
  col = col_fun,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sig_labels[i, j], x, y, gp = gpar(fontsize = 10, col = "grey15", alpha=0.8))
  },
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "ward.D2",
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  heatmap_legend_param = list(title = "Pearson r"),
  row_split = 3
)

# Nice ! Let's remove variables with not much signal
sig_matrix <- p_adj_matrix < 0.05
sig_counts <- colSums(sig_matrix)

# Keep only variables with more than 5 significant correlations
keep_columns <- sig_counts > 5

# Filter matrices
filtered_cor_matrix <- cor_matrix[, keep_columns]
filtered_p_adj_matrix <- p_adj_matrix[, keep_columns]
filtered_sig_labels <- ifelse(filtered_p_adj_matrix > 0.05, "x", "")

Heatmap(
  filtered_cor_matrix,
  col = col_fun,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(filtered_sig_labels[i, j], x, y, gp = gpar(fontsize = 10, col = "grey15", alpha = 0.8))
  },
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "ward.D2",
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  heatmap_legend_param = list(title = "Pearson r"),
  row_split = 3
)

# Add a row annotation for total abundance
total_abund <- colSums(DMSP)
# Log10 transform (add 1 if you want to preserve 0s)
log_total_abund <- log10(total_abund)

# Define breaks in original (linear) space
abund_breaks_raw <- c(10, 100, 1000, 10000, 100000)

# Log-transform them for matching the log10-abundance scale
abund_breaks_log <- log10(abund_breaks_raw)
log_abund_color_fun <- colorRamp2(
  abund_breaks_log,
  viridis(length(abund_breaks_log))
)

# Add legend settings directly to the annotation
row_ha <- rowAnnotation(
  Total_Abundance = log_total_abund,
  col = list(Total_Abundance = log_abund_color_fun),
  annotation_legend_param = list(
    at = abund_breaks_log,
    labels = abund_breaks_raw,  # Show original (non-log) values
    title = "Total Abundance\n(log scale)"
  ),
  width = unit(5, "mm")
)

Heatmap(
  filtered_cor_matrix,
  col = col_fun,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(filtered_sig_labels[i, j], x, y, gp = gpar(fontsize = 10, col = "grey15", alpha = 0.8))
  },
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "ward.D2",
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  heatmap_legend_param = list(title = "Pearson r"),
  left_annotation = row_ha,
  row_split = 3
)

# Adding R2 info
DMSP_RF <- read.table("/Users/emifaure/Documents/ACE/DMSP/DMSP_lyase_RF-extract.tsv", header = TRUE, sep="\t")
DMSP_RF <- DMSP_RF %>% select(AGC_ID, R2_caret)

# Ensure R² values match the filtered_cor_matrix row order
r2_values <- DMSP_RF$R2_caret[match(rownames(filtered_cor_matrix), DMSP_RF$AGC_ID)]
# Choose a gradient for R²
r2_color_fun <- colorRamp2(c(0, 0.25, 0.5, 0.75, 1), plasma(5))

row_ha <- rowAnnotation(
  Total_Abundance = log_total_abund,
  R2 = r2_values,
  col = list(
    Total_Abundance = log_abund_color_fun,
    R2 = r2_color_fun
  ),
  annotation_legend_param = list(
    Total_Abundance = list(
      at = abund_breaks_log,
      labels = abund_breaks_raw,
      title = "Total Abundance"
    ),
    R2 = list(
      at = c(0, 0.5, 1),
      labels = c("0", "0.5", "1"),
      title = "R² from random forest"
    )
  ),
  width = unit(1.0, "cm")
)

ht <- draw(Heatmap(
  filtered_cor_matrix,
  col = col_fun,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(filtered_sig_labels[i, j], x, y, gp = gpar(fontsize = 10, col = "grey15", alpha = 0.8))
  },
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "ward.D2",
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  heatmap_legend_param = list(title = "Pearson r"),
  left_annotation = row_ha,
  row_split = 3
))

# Store the split groups :
split_groups <- row_order(ht)
split_AGC <- lapply(split_groups, function(idx) {
  rownames(filtered_cor_matrix)[idx]
})

# Exploring taxonomy

taxo_dmsp <- fread("Documents/ACE/DMSP/AGC_CONTIG_ORF_EggNOGBestDesc_Kraken-vs-GTDB_DMSPLyase_reformat.tsv", data.table=F, header=F, fill=T, stringsAsFactors=T)
names(taxo_dmsp) = c("AGC_ID","ORF_ID","Contig_ID","EggNOG_Best_Desc","Root", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

summary(taxo_dmsp[which(taxo_dmsp$AGC_ID %in% split_AGC[[1]]),])
summary(taxo_dmsp[which(taxo_dmsp$AGC_ID %in% split_AGC[[2]]),])
summary(taxo_dmsp[which(taxo_dmsp$AGC_ID %in% split_AGC[[3]]),])

# Pie charts :

# Collect genera that are in top 5 of any split
get_top5_per_split <- function(split_ids, taxo_df) {
  taxo_df %>%
    filter(AGC_ID %in% split_ids) %>%
    mutate(Genus = as.character(Genus)) %>%
    filter(Genus != "") %>%
    count(Genus) %>%
    arrange(desc(n)) %>%
    slice_head(n = 5) %>%
    pull(Genus)
}

# Get union of all top 5 genera across splits
top_genera_all_splits <- unique(unlist(lapply(split_AGC, get_top5_per_split, taxo_df = taxo_dmsp)))
top_genera_all_splits
# Taxo :
#root; Bacteria; Pseudomonadota; Alphaproteobacteria; Pelagibacterales; Pelagibacteraceae; Pelagibacter
#root; Bacteria; Pseudomonadota; Alphaproteobacteria; Rhodobacterales; Rhodobacteraceae; Planktomarina
#root; Bacteria; Pseudomonadota; Alphaproteobacteria; Rhodobacterales; Rhodobacteraceae; Amylibacter
#root; Bacteria; Pseudomonadota; Alphaproteobacteria; Rhodobacterales; Rhodobacteraceae; CPC320
#root; Bacteria; Pseudomonadota; Gammaproteobacteria; Pseudomonadales; Nitrincolaceae; ASP10-02a
#root; Bacteria; Pseudomonadota; Alphaproteobacteria; Rhodobacterales; Rhodobacteraceae; GCA-002712045
#root; Bacteria; Pseudomonadota; Alphaproteobacteria; Rhodobacterales; Rhodobacteraceae; MED-G52
#root; Bacteria; Pseudomonadota; Alphaproteobacteria; Rhodobacterales; Rhodobacteraceae; GCA-2697345
#root; Bacteria; Actinomycetota; Acidimicrobiia; Acidimicrobiales; MedAcidi-G1; UBA9410
#root; Bacteria; Pseudomonadota; Gammaproteobacteria; Arenicellales; UBA868; UBA11791
#root; Bacteria; Pseudomonadota; Gammaproteobacteria; Pseudomonadales; Azotimanducaceae_A; UBA9659
#root; Bacteria; Pseudomonadota; Gammaproteobacteria; Arenicellales; UBA868; REDSEA-S09-B13
#root; Bacteria; Pseudomonadota; Alphaproteobacteria; Pelagibacterales; AG-422-B15; CACCCN01

# Generate genus-level colors
genus_colors <- c(
  # Rhodobacteraceae — reds and purples
  "Planktomarina" = "#9e0142",  # dark red
  "Amylibacter" = "#d53e4f",    # red
  "CPC320" = "#d73027",         # orange-red
  "GCA-002712045" = "red",  # burnt orange
  "MED-G52" = "#762a83",         # purple
  "GCA-2697345" = "purple",    # violet
  
  # Pelagibacteraceae
  "Pelagibacter" = "turquoise",    # teal
  
  # Other groups — distinct colors
  "ASP10-02a" = "yellow",       # bright red (Nitrincolaceae)
  "CACCCN01" = "#4575b4",        # blue (AG-422-B15)
  "UBA9410" = "orange",         # aqua (MedAcidi-G1)
  "UBA11791" = "forestgreen",        # light green (UBA868)
  "REDSEA-S09-B13" = "#a6d854",  # green (UBA868)
  "UBA9659" = "black",         # orange (Azotimanducaceae_A)
  
  # Others
  "Others" = "grey50"
)

plot_split_pie <- function(split_ids, taxo_df, genus_colors) {
  
  df <- taxo_df %>% 
    filter(AGC_ID %in% split_ids) %>%
    mutate(Genus = as.character(Genus)) %>%
    count(Genus) %>%
    arrange(desc(n))
  
  # Separate top 5 real genera (excluding empty string)
  top_genera <- df %>% filter(Genus != "") %>% head(5)
  
  # Everything else (incl. unclassified) goes to "Others"
  others <- df %>% filter(!(Genus %in% top_genera$Genus)) %>%
    summarise(n = sum(n)) %>%
    pull(n)
  
  pie_df <- rbind(top_genera, data.frame(Genus = "Others", n = others))
  
  pie_df$Fraction <- pie_df$n / sum(pie_df$n)
  pie_df$Label <- paste0(pie_df$Genus, " (", round(100 * pie_df$Fraction), "%)")
  
  ggplot(pie_df, aes(x = "", y = Fraction, fill = Genus)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(legend.title = element_blank()) +
    scale_fill_manual(values = genus_colors)
}

plot_split_pie(split_AGC[[1]], taxo_dmsp, genus_colors)
plot_split_pie(split_AGC[[2]], taxo_dmsp, genus_colors)
plot_split_pie(split_AGC[[3]], taxo_dmsp, genus_colors)

# Get the legend :
barplot(rep(1, length(genus_colors)),
        col = genus_colors,
        names.arg = names(genus_colors),
        las = 2)

########## Trying Redundancy analysis :

DMSP_hel <- decostand(DMSP,method = "hellinger")

RDA <- rda(DMSP_hel ~ ., data=meta_DMSP)
#summary(RDA)
RsquareAdj(RDA) #81.4%
anova.cca(RDA) #Significant
anova.cca(RDA, by="axis") # 5 significant

# We need to select environmental variables before going further with plots
set.seed(1994)
rda.step.both = ordistep(rda(DMSP_hel ~1,data=meta_DMSP),scope=formula(RDA),direction="both",steps=100)
res.rda = rda(formula=formula(rda.step.both), data=meta_DMSP)
RsquareAdj(res.rda) #80.2
#Test significance
# Full model
anova.rda <- anova(res.rda) # Significant 0.001
# By axis
anova.rda.ax <- anova(res.rda, by = "axis", permutations = how(nperm = 499))
# 5 significant

scores=scoresRDA(res.rda)

samples = scores$Sites
samples = as.data.frame(samples)
samples = merge(samples,meta_DMSP, by="row.names")
row.names(samples)=samples[,1]
samples=samples[,-1]
# Add group variable to color points in the right colorway :
Group_Var=Group_Var[Group_Var$ACE_seq_name %in% row.names(samples),]
row.names(Group_Var)=Group_Var$ACE_seq_name
Group_Var=Group_Var[,-1]
samples = merge(samples,Group_Var,by="row.names")
rownames(samples) = samples[,1]
samples = samples[,-1]

# Retrieve AGC scores
AGCScore <- scores$Species
AGCScore = as.data.frame(AGCScore)

enviscore=scores$Biplot

#Plot the triplot
ggplot() +   geom_hline(yintercept = 0, linetype='dotted') +
  geom_vline(xintercept = 0, linetype='dotted') +
  labs(x = paste0("RDA1 (",
                  round(100 * scores$Eigenval.rel[1], 1), "%)"),
       y = paste0("RDA2 (",
                  round(100 * scores$Eigenval.rel[2], 1), "%)")) +
  theme(plot.title=element_text(hjust=0.5)) +
  geom_point(aes(samples[,1], samples[,2], col=samples$Water_mass_simplified.x, shape=samples$MertzGlacier.x), size = 2) +
  geom_segment(data=as.data.frame(enviscore),aes(xend = RDA1, yend = RDA2),x=0,y=0,linewidth = 0.5, linetype="F1",color = 'steelblue4',arrow = arrow(length = unit(0.2,"cm"))) +
  geom_text(aes(enviscore[,1]*1.05, enviscore[,2]*1.05, label=rownames(enviscore)), color="steelblue4") +
  geom_segment(data=as.data.frame(AGCScore),aes(xend = RDA1, yend = RDA2),x=0,y=0,linewidth = 0.5,color = 'indianred2',arrow = arrow(length = unit(0.2,"cm"))) +
  geom_text(data=as.data.frame(AGCScore[which(abs(AGCScore$RDA1)>0.5 | abs(AGCScore$RDA2)>0.5),]),aes(x = RDA1, y = RDA2, label=rownames(AGCScore[which(abs(AGCScore$RDA1)>0.5 | abs(AGCScore$RDA2)>0.5),])), color='indianred2')+
  scale_colour_manual(values = levels(samples$Water_mass_simplified.col)) +
  theme_bw() 
