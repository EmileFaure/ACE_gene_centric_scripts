library(tidyverse)
library(ggVennDiagram)
library(viridis)
library(ComplexHeatmap)
library(circlize)
set.seed(1994)

# AGCs
Functions_22=read.table("/Users/emifaure/Documents/Gradient_Explo/Post_Review/AGCCountsPerEggNOG_GTDBandUNIREF_SAR11/CAG_22_AGCCountsPerEggNOG_GTDBandUNIREF_SAR11",sep="\t", header=F, stringsAsFactors=T, fill=T, quote="")
names(Functions_22)=c("GeneCount","Function")
Functions_22$CAG="CAG_22"
Functions_35=read.table("/Users/emifaure/Documents/Gradient_Explo/Post_Review/AGCCountsPerEggNOG_GTDBandUNIREF_SAR11/CAG_35_AGCCountsPerEggNOG_GTDBandUNIREF_SAR11",sep="\t", header=F, stringsAsFactors=T, fill=T, quote="")
names(Functions_35)=c("GeneCount","Function")
Functions_35$CAG="CAG_35"
Functions_83=read.table("/Users/emifaure/Documents/Gradient_Explo/Post_Review/AGCCountsPerEggNOG_GTDBandUNIREF_SAR11/CAG_83_AGCCountsPerEggNOG_GTDBandUNIREF_SAR11",sep="\t", header=F, stringsAsFactors=T, fill=T, quote="")
names(Functions_83)=c("GeneCount","Function")
Functions_83$CAG="CAG_83"

Functions=rbind(Functions_22,Functions_35,Functions_83)

length(which(table(Functions$Function)==3)) # 1061
length(which(table(Functions$Function)==2)) # 384
length(which(table(Functions$Function)==1)) # 600

length(unique(Functions$Function)) # 2045

Unique_Func=Functions[which(Functions$Function %in% names(which(table(Functions$Function)==1))),]
summary(Unique_Func)

Unique_Func[Unique_Func$CAG=="CAG_83" & Unique_Func$GeneCount >= 5,c(1,2)]
Unique_Func[Unique_Func$CAG=="CAG_35" & Unique_Func$GeneCount > 5,c(1,2)]
Unique_Func[Unique_Func$CAG=="CAG_22"& Unique_Func$GeneCount >= 5,c(1,2)]

Data_Venn=split(as.character(Functions$Function), f = Functions$CAG)
Venn_EggNog=ggVennDiagram(Data_Venn) #+ scale_fill_viridis()
Venn_EggNog

# Let's try a different representation
Functions_22_HM=Functions_22[,c(1,2)]
names(Functions_22_HM)=c("GeneCount_CAG22","Function")
Functions_35_HM=Functions_35[,c(1,2)]
names(Functions_35_HM)=c("GeneCount_CAG35","Function")
Functions_heatmap=merge(Functions_22_HM,Functions_35_HM,all = T)
Functions_83_HM=Functions_83[,c(1,2)]
names(Functions_83_HM)=c("GeneCount_CAG83","Function")
Functions_heatmap=merge(Functions_heatmap,Functions_83_HM,all = T)
Functions_heatmap[is.na(Functions_heatmap)]<-0
row.names(Functions_heatmap)=Functions_heatmap$Function
Functions_heatmap=Functions_heatmap[,-1]

# Remove the non-annotated :
Functions_heatmap[1:3,1:3]
Functions_heatmap=Functions_heatmap[-c(1:3),]
# Transpose and remove names for graphic readability :
Functions_heatmap_transposed=t(Functions_heatmap)
# To be launched for better graphical rendering :
colnames(Functions_heatmap_transposed)=NULL
row.names(Functions_heatmap_transposed)=c("CAG22","CAG35","CAG83")

# Log rep :
col_fun = colorRamp2(c(0,2.5,5), c("grey20", "indianred3", "gold"))
Hm=Heatmap(log(Functions_heatmap_transposed+1),
           name="Scaled normalized\nSAR11-related AGC counts\nper function",
           row_title="CAGs of interest",
           row_order=c(1,2,3),
           #column_title=c("Increased presence in CAG35\n(Ubiquitous, positively correlated\nwith temperature)","Shared presence in\nCAG83 and CAG35","Increased presence in CAG 83\n(Subtropical waters)","Shared presence in\nCAG83 and CAG22","Increased presence in CAG22\n(Antarctic surface waters)"),
           column_dend_height = unit(2, "cm"),
           clustering_distance_columns = "euclidean",
           clustering_method_columns = "ward.D2",
           column_split = 6,
           column_title_gp = gpar(fontsize=10, fontface="italic"),
           column_names_gp = gpar(fontsize=4),
           column_title_side = "bottom",
           column_title_rot = 90,
           col=col_fun)
Hm=draw(Hm)

col_fun = colorRamp2(c(0,0.5,1), c("grey20", "indianred3", "gold"))
Hm=Heatmap(sweep(Functions_heatmap_transposed, 2, colSums(Functions_heatmap_transposed), FUN = "/"),
           name="Proportion of AGC\nin each CAG\nper function\n(only SAR11-related AGCs)",
           row_title="CAGs of interest",
           row_order=c(1,2,3),
           column_title=c("Specific CAG35\n(Ubiquitous, positively correlated\nwith temperature)","Specific CAG83\n(Subtropical waters)","Increased presence in CAG 35", "Shared presence in \nCAG35 and CAG22", "Homogeneous presence", "Shared presence in\nCAG83 and CAG35","Specific to CAG 22\n(Antarctic surface waters)", "Increased presence in CAG22"),
           column_dend_height = unit(2, "cm"),
           clustering_distance_columns = "euclidean",
           clustering_method_columns = "ward.D2",
           column_split = 8,
           column_title_gp = gpar(fontsize=10, fontface="italic"),
           column_names_gp = gpar(fontsize=4),
           column_title_side = "bottom",
           column_title_rot = 90,
           col=col_fun)
Hm=draw(Hm)

#### Exploration of AGCs annotated to functions specific of the different CAGs

########## Load metadata ############

meta <- readRDS("/Users/emifaure/Documents/ACE/Metadata_GeneMat_KNN_MertzTF_1099Corr.rds")
Group_Var=readRDS("/Users/emifaure/Documents/ACE/metadata_metaG_grouping_variables.rds")

########## Load AGC abundance ############
focalfunc <- "NikM"
Enzyme <- read.table(paste0("/Users/emifaure/Documents/Gradient_Explo/GeneMats_AGC_SAR/GM_AGN_TMAX60GQ_transposed_NZVuniquecut20_FreeLiving_",focalfunc,".tsv"), header = TRUE, sep="\t")

# Select metadata from samples :
meta_Enzyme <- meta[which(meta$ACE_seq_name %in% row.names(Enzyme)),]
rownames(meta_Enzyme) <- meta_Enzyme$ACE_seq_name
meta_Enzyme <- meta_Enzyme[,-1]
# Match row order
meta_Enzyme <- meta_Enzyme[match(row.names(Enzyme), row.names(meta_Enzyme)),]
#extract numeric for correlation computations :
meta_Enzyme_quanti <- meta_Enzyme %>% select(where(is.numeric))
#The sign of some latent variables needs to be inversed to be interpreted as in the Landwehr et al. study
#Mail from Landwehr : "We have chosen to change the signs of some of the LVs that where produced by the algorithm, in order to make them align with certain physical processes"
#Also note that here LV are from 0 to 13, like in the data provided by Landwehr et al., need to add one to each LV number to match their study.
meta_Enzyme_quanti$LV3_mean = -meta_Enzyme_quanti$LV3_mean
meta_Enzyme_quanti$LV12_mean = -meta_Enzyme_quanti$LV12_mean
meta_Enzyme_quanti$LV4_mean = -meta_Enzyme_quanti$LV4_mean
meta_Enzyme_quanti$LV11_mean = -meta_Enzyme_quanti$LV11_mean
meta_Enzyme_quanti$LV10_mean = -meta_Enzyme_quanti$LV10_mean
meta_Enzyme_quanti$LV5_mean = -meta_Enzyme_quanti$LV5_mean
meta_Enzyme_quanti$LV9_mean = -meta_Enzyme_quanti$LV9_mean

# Initialize matrices
n_abund <- ncol(Enzyme)
n_env <- ncol(meta_Enzyme_quanti)

cor_matrix <- matrix(NA, nrow = n_abund, ncol = n_env)
p_matrix <- matrix(NA, nrow = n_abund, ncol = n_env)

# Loop over all combinations of abundance vs environment
for (i in 1:n_abund) {
  for (j in 1:n_env) {
    test <- cor.test(Enzyme[, i], meta_Enzyme_quanti[, j], method = "pearson")
    cor_matrix[i, j] <- test$estimate
    p_matrix[i, j] <- test$p.value
  }
}

# Assign row/column names for clarity
rownames(cor_matrix) <- colnames(Enzyme)
colnames(cor_matrix) <- colnames(meta_Enzyme_quanti)

rownames(p_matrix) <- colnames(Enzyme)
colnames(p_matrix) <- colnames(meta_Enzyme_quanti)

# Adjust
p_adj_vector <- p.adjust(as.vector(p_matrix), method = "BH")
p_adj_matrix <- matrix(p_adj_vector, nrow = n_abund, ncol = n_env)
rownames(p_adj_matrix) <- colnames(Enzyme)
colnames(p_adj_matrix) <- colnames(meta_Enzyme_quanti)

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
  heatmap_legend_param = list(title = "Pearson r")
)

# Nice ! Let's remove variables with not much signal
sig_matrix <- p_adj_matrix < 0.05
sig_counts <- colSums(sig_matrix)

# Keep only variables with more than 5 significant correlations
keep_columns <- sig_counts > 4

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
  heatmap_legend_param = list(title = "Pearson r")
)

# Add a row annotation for total abundance
total_abund <- colSums(Enzyme)
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
  left_annotation = row_ha
)


# Adding % of SAR11 info
Enzyme_SAR11 <- read.table(paste0("/Users/emifaure/Documents/Gradient_Explo/GeneMats_AGC_SAR/AGC_",focalfunc,".stats"), header = FALSE, sep="\t")
names(Enzyme_SAR11) <- c("AGC_ID","AGC_Size","SAR11_ORFs")
Enzyme_SAR11$SAR11prop <- Enzyme_SAR11$SAR11_ORFs / Enzyme_SAR11$AGC_Size
Enzyme_SAR11 <- Enzyme_SAR11 %>% select(AGC_ID, SAR11prop)

# Ensure R² values match the filtered_cor_matrix row order
SAR11_values <- Enzyme_SAR11$SAR11prop[match(rownames(filtered_cor_matrix), Enzyme_SAR11$AGC_ID)]
# Choose a gradient for R²
SAR11_color_fun <- colorRamp2(c(0, 0.25, 0.5, 0.75, 1), plasma(5))

row_ha <- rowAnnotation(
  Total_Abundance = log_total_abund,
  SAR11_prop = SAR11_values,
  col = list(
    Total_Abundance = log_abund_color_fun,
    SAR11_prop = SAR11_color_fun
  ),
  annotation_legend_param = list(
    Total_Abundance = list(
      at = abund_breaks_log,
      labels = abund_breaks_raw,
      title = "Total Abundance"
    ),
    SAR11_prop = list(
      at = c(0, 0.5, 1),
      labels = c("0", "0.5", "1"),
      title = "Prop. of SAR11 ORFs"
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
  left_annotation = row_ha
))

# Building one heatmap per CAG
  # CAG22 :
focalfunc <- "NikM"
NikM <- read.table(paste0("/Users/emifaure/Documents/Gradient_Explo/GeneMats_AGC_SAR/GM_AGN_TMAX60GQ_transposed_NZVuniquecut20_FreeLiving_",focalfunc,".tsv"), header = TRUE, sep="\t")
focalfunc <- "PNDOxido"
PNDOxido <- read.table(paste0("/Users/emifaure/Documents/Gradient_Explo/GeneMats_AGC_SAR/GM_AGN_TMAX60GQ_transposed_NZVuniquecut20_FreeLiving_",focalfunc,".tsv"), header = TRUE, sep="\t")
focalfunc <- "LysM"
LysM <- read.table(paste0("/Users/emifaure/Documents/Gradient_Explo/GeneMats_AGC_SAR/GM_AGN_TMAX60GQ_transposed_NZVuniquecut20_FreeLiving_",focalfunc,".tsv"), header = TRUE, sep="\t")

Enzyme <- merge(NikM, PNDOxido, by="row.names")
row.names(Enzyme) <- Enzyme[,1]
Enzyme <- Enzyme[,-1]
Enzyme <- merge(Enzyme, LysM, by="row.names")
row.names(Enzyme) <- Enzyme[,1]
Enzyme <- Enzyme[,-1]

# Select metadata from samples :
meta_Enzyme <- meta[which(meta$ACE_seq_name %in% row.names(Enzyme)),]
rownames(meta_Enzyme) <- meta_Enzyme$ACE_seq_name
meta_Enzyme <- meta_Enzyme[,-1]
# Match row order
meta_Enzyme <- meta_Enzyme[match(row.names(Enzyme), row.names(meta_Enzyme)),]
#extract numeric for correlation computations :
meta_Enzyme_quanti <- meta_Enzyme %>% select(where(is.numeric))
#The sign of some latent variables needs to be inversed to be interpreted as in the Landwehr et al. study
#Mail from Landwehr : "We have chosen to change the signs of some of the LVs that where produced by the algorithm, in order to make them align with certain physical processes"
#Also note that here LV are from 0 to 13, like in the data provided by Landwehr et al., need to add one to each LV number to match their study.
meta_Enzyme_quanti$LV3_mean = -meta_Enzyme_quanti$LV3_mean
meta_Enzyme_quanti$LV12_mean = -meta_Enzyme_quanti$LV12_mean
meta_Enzyme_quanti$LV4_mean = -meta_Enzyme_quanti$LV4_mean
meta_Enzyme_quanti$LV11_mean = -meta_Enzyme_quanti$LV11_mean
meta_Enzyme_quanti$LV10_mean = -meta_Enzyme_quanti$LV10_mean
meta_Enzyme_quanti$LV5_mean = -meta_Enzyme_quanti$LV5_mean
meta_Enzyme_quanti$LV9_mean = -meta_Enzyme_quanti$LV9_mean

# Initialize matrices
n_abund <- ncol(Enzyme)
n_env <- ncol(meta_Enzyme_quanti)

cor_matrix <- matrix(NA, nrow = n_abund, ncol = n_env)
p_matrix <- matrix(NA, nrow = n_abund, ncol = n_env)

# Loop over all combinations of abundance vs environment
for (i in 1:n_abund) {
  for (j in 1:n_env) {
    test <- cor.test(Enzyme[, i], meta_Enzyme_quanti[, j], method = "pearson")
    cor_matrix[i, j] <- test$estimate
    p_matrix[i, j] <- test$p.value
  }
}

# Assign row/column names for clarity
rownames(cor_matrix) <- colnames(Enzyme)
colnames(cor_matrix) <- colnames(meta_Enzyme_quanti)

rownames(p_matrix) <- colnames(Enzyme)
colnames(p_matrix) <- colnames(meta_Enzyme_quanti)

# Adjust
p_adj_vector <- p.adjust(as.vector(p_matrix), method = "BH")
p_adj_matrix <- matrix(p_adj_vector, nrow = n_abund, ncol = n_env)
rownames(p_adj_matrix) <- colnames(Enzyme)
colnames(p_adj_matrix) <- colnames(meta_Enzyme_quanti)

# Define color scale
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Create matrix of significance labels (e.g., "*" for p < 0.05)
sig_labels <- ifelse(p_adj_matrix > 0.05, "x", "")

# Split by function
NikM_cols <- colnames(NikM)
PNDOxido_cols <- colnames(PNDOxido)
LysM_cols <- colnames(LysM)
row_group <- c(
  rep("NikM", length(NikM_cols)),
  rep("PNDOxido", length(PNDOxido_cols)),
  rep("LysM", length(LysM_cols))
)

# Let's remove variables with not much signal
sig_matrix <- p_adj_matrix < 0.05
sig_counts <- colSums(sig_matrix)

# Keep only variables with more than 10 significant correlations
keep_columns <- sig_counts > 10

# Filter matrices
filtered_cor_matrix <- cor_matrix[, keep_columns]
filtered_p_adj_matrix <- p_adj_matrix[, keep_columns]
filtered_sig_labels <- ifelse(filtered_p_adj_matrix > 0.05, "x", "")

# Add a row annotation for total abundance
total_abund <- colSums(Enzyme)
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

# Adding % of SAR11 info
focalfunc <- "NikM"
NikM_SAR11 <- read.table(paste0("/Users/emifaure/Documents/Gradient_Explo/GeneMats_AGC_SAR/AGC_",focalfunc,".stats"), header = FALSE, sep="\t")
focalfunc <- "PNDOxido"
PNDOxido_SAR11 <- read.table(paste0("/Users/emifaure/Documents/Gradient_Explo/GeneMats_AGC_SAR/AGC_",focalfunc,".stats"), header = FALSE, sep="\t")
focalfunc <- "LysM"
LysM_SAR11 <- read.table(paste0("/Users/emifaure/Documents/Gradient_Explo/GeneMats_AGC_SAR/AGC_",focalfunc,".stats"), header = FALSE, sep="\t")
Enzyme_SAR11 <- as.data.frame(rbind(NikM_SAR11, PNDOxido_SAR11, LysM_SAR11))
names(Enzyme_SAR11) <- c("AGC_ID","AGC_Size","SAR11_ORFs")
Enzyme_SAR11$SAR11prop <- Enzyme_SAR11$SAR11_ORFs / Enzyme_SAR11$AGC_Size
Enzyme_SAR11 <- Enzyme_SAR11 %>% select(AGC_ID, SAR11prop)

# Ensure R² values match the filtered_cor_matrix row order
SAR11_values <- Enzyme_SAR11$SAR11prop[match(rownames(filtered_cor_matrix), Enzyme_SAR11$AGC_ID)]
# Choose a gradient for R²
SAR11_color_fun <- colorRamp2(c(0, 0.25, 0.5, 0.75, 1), plasma(5))

row_ha <- rowAnnotation(
  Total_Abundance = log_total_abund,
  SAR11_prop = SAR11_values,
  col = list(
    Total_Abundance = log_abund_color_fun,
    SAR11_prop = SAR11_color_fun
  ),
  annotation_legend_param = list(
    Total_Abundance = list(
      at = abund_breaks_log,
      labels = abund_breaks_raw,
      title = "Total Abundance"
    ),
    SAR11_prop = list(
      at = c(0, 0.5, 1),
      labels = c("0", "0.5", "1"),
      title = "Prop. of SAR11 ORFs"
    )
  ),
  width = unit(1.0, "cm")
)

# Further select and organize variables :
selectedvars <-  c("Nitrite_µmol.L","dO2_µmol.kg","PON.µmol","POC_µM","Bsi_µM","dCd_nmol.kg","dZn_nmol.kg",          
                   "Depth_m","d13C","dFe_nmol.kg","Salinity_PSS78")
filtered_cor_matrix <- filtered_cor_matrix[, selectedvars]
filtered_p_adj_matrix <- filtered_p_adj_matrix[, selectedvars]
filtered_sig_labels <- ifelse(filtered_p_adj_matrix > 0.05, "x", "")

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
  row_split = row_group,
  cluster_row_slices = FALSE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  heatmap_legend_param = list(title = "Pearson r"),
  left_annotation = row_ha
))

###################################

# CAG83 :
focalfunc <- "ABCdipep"
ABCdipep <- read.table(paste0("/Users/emifaure/Documents/Gradient_Explo/GeneMats_AGC_SAR/GM_AGN_TMAX60GQ_transposed_NZVuniquecut20_FreeLiving_",focalfunc,".tsv"), header = TRUE, sep="\t")
focalfunc <- "NiABC"
NiABC <- read.table(paste0("/Users/emifaure/Documents/Gradient_Explo/GeneMats_AGC_SAR/GM_AGN_TMAX60GQ_transposed_NZVuniquecut20_FreeLiving_",focalfunc,".tsv"), header = TRUE, sep="\t")
focalfunc <- "Zn"
Zn <- read.table(paste0("/Users/emifaure/Documents/Gradient_Explo/GeneMats_AGC_SAR/GM_AGN_TMAX60GQ_transposed_NZVuniquecut20_FreeLiving_",focalfunc,".tsv"), header = TRUE, sep="\t")
focalfunc <- "GlycineBetaine"
GlycineBetaine <- read.table(paste0("/Users/emifaure/Documents/Gradient_Explo/GeneMats_AGC_SAR/GM_AGN_TMAX60GQ_transposed_NZVuniquecut20_FreeLiving_",focalfunc,".tsv"), header = TRUE, sep="\t")

Enzyme <- merge(ABCdipep, NiABC, by="row.names")
row.names(Enzyme) <- Enzyme[,1]
Enzyme <- Enzyme[,-1]
Enzyme <- merge(Enzyme, Zn, by="row.names")
row.names(Enzyme) <- Enzyme[,1]
Enzyme <- Enzyme[,-1]
Enzyme <- merge(Enzyme, GlycineBetaine, by="row.names")
row.names(Enzyme) <- Enzyme[,1]
Enzyme <- Enzyme[,-1]

# Select metadata from samples :
meta_Enzyme <- meta[which(meta$ACE_seq_name %in% row.names(Enzyme)),]
rownames(meta_Enzyme) <- meta_Enzyme$ACE_seq_name
meta_Enzyme <- meta_Enzyme[,-1]
# Match row order
meta_Enzyme <- meta_Enzyme[match(row.names(Enzyme), row.names(meta_Enzyme)),]
#extract numeric for correlation computations :
meta_Enzyme_quanti <- meta_Enzyme %>% select(where(is.numeric))
#The sign of some latent variables needs to be inversed to be interpreted as in the Landwehr et al. study
#Mail from Landwehr : "We have chosen to change the signs of some of the LVs that where produced by the algorithm, in order to make them align with certain physical processes"
#Also note that here LV are from 0 to 13, like in the data provided by Landwehr et al., need to add one to each LV number to match their study.
meta_Enzyme_quanti$LV3_mean = -meta_Enzyme_quanti$LV3_mean
meta_Enzyme_quanti$LV12_mean = -meta_Enzyme_quanti$LV12_mean
meta_Enzyme_quanti$LV4_mean = -meta_Enzyme_quanti$LV4_mean
meta_Enzyme_quanti$LV11_mean = -meta_Enzyme_quanti$LV11_mean
meta_Enzyme_quanti$LV10_mean = -meta_Enzyme_quanti$LV10_mean
meta_Enzyme_quanti$LV5_mean = -meta_Enzyme_quanti$LV5_mean
meta_Enzyme_quanti$LV9_mean = -meta_Enzyme_quanti$LV9_mean

# Initialize matrices
n_abund <- ncol(Enzyme)
n_env <- ncol(meta_Enzyme_quanti)

cor_matrix <- matrix(NA, nrow = n_abund, ncol = n_env)
p_matrix <- matrix(NA, nrow = n_abund, ncol = n_env)

# Loop over all combinations of abundance vs environment
for (i in 1:n_abund) {
  for (j in 1:n_env) {
    test <- cor.test(Enzyme[, i], meta_Enzyme_quanti[, j], method = "pearson")
    cor_matrix[i, j] <- test$estimate
    p_matrix[i, j] <- test$p.value
  }
}

# Assign row/column names for clarity
rownames(cor_matrix) <- colnames(Enzyme)
colnames(cor_matrix) <- colnames(meta_Enzyme_quanti)

rownames(p_matrix) <- colnames(Enzyme)
colnames(p_matrix) <- colnames(meta_Enzyme_quanti)

# Adjust
p_adj_vector <- p.adjust(as.vector(p_matrix), method = "BH")
p_adj_matrix <- matrix(p_adj_vector, nrow = n_abund, ncol = n_env)
rownames(p_adj_matrix) <- colnames(Enzyme)
colnames(p_adj_matrix) <- colnames(meta_Enzyme_quanti)

# Define color scale
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Create matrix of significance labels (e.g., "*" for p < 0.05)
sig_labels <- ifelse(p_adj_matrix > 0.05, "x", "")

# Split by function
ABCdipep_cols <- colnames(ABCdipep)
NiABC_cols <- colnames(NiABC)
Zn_cols <- colnames(Zn)
GlycineBetaine_cols <- colnames(GlycineBetaine)
row_group <- c(
  rep("ABCdipep", length(ABCdipep_cols)),
  rep("NiABC", length(NiABC_cols)),
  rep("Zn", length(Zn_cols)),
  rep("GlycineBetaine", length(GlycineBetaine_cols))
)

# Let's remove variables with not much signal
sig_matrix <- p_adj_matrix < 0.05
sig_counts <- colSums(sig_matrix)

# Keep only variables with more than 10 significant correlations
keep_columns <- sig_counts > 10

# Filter matrices
filtered_cor_matrix <- cor_matrix[, keep_columns]
filtered_p_adj_matrix <- p_adj_matrix[, keep_columns]
filtered_sig_labels <- ifelse(filtered_p_adj_matrix > 0.05, "x", "")

# Add a row annotation for total abundance
total_abund <- colSums(Enzyme)
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

# Adding % of SAR11 info
focalfunc <- "ABCdipep"
ABCdipep_SAR11 <- read.table(paste0("/Users/emifaure/Documents/Gradient_Explo/GeneMats_AGC_SAR/AGC_",focalfunc,".stats"), header = FALSE, sep="\t")
focalfunc <- "NiABC"
NiABC_SAR11 <- read.table(paste0("/Users/emifaure/Documents/Gradient_Explo/GeneMats_AGC_SAR/AGC_",focalfunc,".stats"), header = FALSE, sep="\t")
focalfunc <- "Zn"
Zn_SAR11 <- read.table(paste0("/Users/emifaure/Documents/Gradient_Explo/GeneMats_AGC_SAR/AGC_",focalfunc,".stats"), header = FALSE, sep="\t")
focalfunc <- "GlycineBetaine"
GlycineBetaine_SAR11 <- read.table(paste0("/Users/emifaure/Documents/Gradient_Explo/GeneMats_AGC_SAR/AGC_",focalfunc,".stats"), header = FALSE, sep="\t")
Enzyme_SAR11 <- as.data.frame(rbind(ABCdipep_SAR11, NiABC_SAR11, Zn_SAR11, GlycineBetaine_SAR11))
names(Enzyme_SAR11) <- c("AGC_ID","AGC_Size","SAR11_ORFs")
Enzyme_SAR11$SAR11prop <- Enzyme_SAR11$SAR11_ORFs / Enzyme_SAR11$AGC_Size
Enzyme_SAR11 <- Enzyme_SAR11 %>% select(AGC_ID, SAR11prop)

# Ensure R² values match the filtered_cor_matrix row order
SAR11_values <- Enzyme_SAR11$SAR11prop[match(rownames(filtered_cor_matrix), Enzyme_SAR11$AGC_ID)]
# Choose a gradient for R²
SAR11_color_fun <- colorRamp2(c(0, 0.25, 0.5, 0.75, 1), plasma(5))

row_ha <- rowAnnotation(
  Total_Abundance = log_total_abund,
  SAR11_prop = SAR11_values,
  col = list(
    Total_Abundance = log_abund_color_fun,
    SAR11_prop = SAR11_color_fun
  ),
  annotation_legend_param = list(
    Total_Abundance = list(
      at = abund_breaks_log,
      labels = abund_breaks_raw,
      title = "Total Abundance"
    ),
    SAR11_prop = list(
      at = c(0, 0.5, 1),
      labels = c("0", "0.5", "1"),
      title = "Prop. of SAR11 ORFs"
    )
  ),
  width = unit(1.0, "cm")
)

# Further select and organize variables :
selectedvars <-  c("Nitrate_µmol.L","Phosphate_µmol.L","dCd_nmol.kg","dZn_nmol.kg","dO2_µmol.kg",          
                  "Salinity_PSS78","Latitude","d13C","Temperature_ITS90")
filtered_cor_matrix <- filtered_cor_matrix[, selectedvars]
filtered_p_adj_matrix <- filtered_p_adj_matrix[, selectedvars]
filtered_sig_labels <- ifelse(filtered_p_adj_matrix > 0.05, "x", "")

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
  row_split = row_group,
  cluster_row_slices = FALSE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  heatmap_legend_param = list(title = "Pearson r"),
  left_annotation = row_ha
))

############### THRASH  #######################


# ORFs
Functions_22=read.table("/Users/emifaure/Documents/Gradient_Explo/Post_Review/ORFCountsPerEggNOG_GTDBandUNIREF_SAR11/CAG_22_ORFCountsPerEggNOG_GTDBandUNIREF_SAR11",sep="\t", header=F, stringsAsFactors=T, fill=T, quote="")
names(Functions_22)=c("GeneCount","Function")
Functions_22$CAG="CAG_22"
Functions_35=read.table("/Users/emifaure/Documents/Gradient_Explo/Post_Review/ORFCountsPerEggNOG_GTDBandUNIREF_SAR11/CAG_35_ORFCountsPerEggNOG_GTDBandUNIREF_SAR11",sep="\t", header=F, stringsAsFactors=T, fill=T, quote="")
names(Functions_35)=c("GeneCount","Function")
Functions_35$CAG="CAG_35"
Functions_83=read.table("/Users/emifaure/Documents/Gradient_Explo/Post_Review/ORFCountsPerEggNOG_GTDBandUNIREF_SAR11/CAG_83_ORFCountsPerEggNOG_GTDBandUNIREF_SAR11",sep="\t", header=F, stringsAsFactors=T, fill=T, quote="")
names(Functions_83)=c("GeneCount","Function")
Functions_83$CAG="CAG_83"

Functions=rbind(Functions_22,Functions_35,Functions_83)

length(which(table(Functions$Function)==3)) # 1061
length(which(table(Functions$Function)==2)) # 384
length(which(table(Functions$Function)==1)) # 600

length(unique(Functions$Function)) # 2045

Unique_Func=Functions[which(Functions$Function %in% names(which(table(Functions$Function)==1))),]
summary(Unique_Func)

Unique_Func[Unique_Func$CAG=="CAG_83" & Unique_Func$GeneCount >= 100,c(1,2)]
Unique_Func[Unique_Func$CAG=="CAG_35" & Unique_Func$GeneCount > 100,c(1,2)]
Unique_Func[Unique_Func$CAG=="CAG_22"& Unique_Func$GeneCount >= 100,c(1,2)]

Data_Venn=split(as.character(Functions$Function), f = Functions$CAG)
Venn_EggNog=ggVennDiagram(Data_Venn) + scale_fill_viridis()
Venn_EggNog
