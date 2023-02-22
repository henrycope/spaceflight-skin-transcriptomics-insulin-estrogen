# Load libraries
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(tidyr)

# Set output directory
output_dir <- "plots"

# Create output directory if it does not exist
if (!file.exists(output_dir)) 
  dir.create(output_dir)

# Load FGSEA results from file
df_gsea <- read.csv("fgsea.csv", header = T)

# Create normalized enrichment score matrix
simple_df <- df_gsea[, c("pathway", "NES", "subset")]
NES_df <- pivot_wider(simple_df, names_from = pathway, values_from = NES)
NES_matrix <- data.matrix(NES_df[, 2:ncol(NES_df)])
rownames(NES_matrix) <- pull(NES_df[1], subset)

# Load metadata
meta_data_df <- read.csv("subsets.csv", row.names = 1, encoding = "UTF-8")

# Heatmap annotation bar colours
anno_colors <- list(
  Mission = c("MHU-2" = "#1D91C0", "RR-5" = "#F7FCFD", "RR-7" = "#CB181D"),
  Strain = c("C57BL" = "#0072b2", "BALB" = "#914770", "C3H" = "#7f7f7f"),
  Diet = c("JC" = "#762A83", "JCwFOS" = "#FFEE99", "NuRFB" = "#225555"),
  Euthanasia = c("LAR" = "#EE8866", "non-LAR" = "#EEDD88"),
  Duration = c("25" = "#f419ec", "30" = "#4419f4", "75" = "#19f4f2"),
  Tissue = c("Dorsal" = "#222255", "Femoral" = "#663333"),
  Age = c("9" = "#44BB99", "11" = "#99DDFF", "30" = "#77AADD"),
  Recovery = c("1" = "#F26E01", "30" = "#E34A27", "0" = "#FFFFFF")
)

# Heatmap annotation bar
leftha <- rowAnnotation(
  Mission = meta_data_df$mission,
  Strain = meta_data_df$strain,
  Diet = meta_data_df$diet,
  Duration = meta_data_df$duration,
  Age = meta_data_df$age,
  Tissue = meta_data_df$tissue,
  Euthanasia = meta_data_df$euthanasia,
  Recovery = meta_data_df$recovery,
  col = anno_colors,
  annotation_legend_param = list(Age = list(at = c(9, 11, 30)), Duration = list(at = c(25, 30, 75))),
  gap = unit(0, "points"),
  border = T,
  annotation_name_gp = gpar(cex = 1, col = "blue"),
  width = unit(4, "cm"),
  simple_anno_size_adjust = TRUE
)

# Heatmap colour scale
col_fun <- colorRamp2(c(min(NES_matrix, na.rm = T), 0, max(NES_matrix, na.rm = T)), c("orange", "white", "green"))

# Function to add levels of significance
cell_fun <- function(j, i, x, y, width, height, fill) {
  padj_value <- df_gsea %>%
    filter(subset == rownames(NES_matrix)[i],
           pathway == colnames(NES_matrix)[j]) %>%
    pull(padj)
  
  if (is.na(padj_value))
    return()
  
  if (padj_value <= 0.01)
    grid.text("***", x, y, gp = gpar(fontsize = 16))
  else if (padj_value <= 0.05)
    grid.text("**", x, y, gp = gpar(fontsize = 16))
  else if (padj_value <= 0.1)
    grid.text("*", x, y, gp = gpar(fontsize = 16))
}

# Legend for significance levels
significance_lgd <- Legend(
  labels = c("≤ 0.10", "≤ 0.05", "≤ 0.01"),
  title = "Significance (FDR)",
  type = "points",
  pch = c("*", "**", "***"),
  background = "white"
)

# Build heatmap
NES_heatmap <- Heatmap(
  NES_matrix,
  name = "NES",
  col = col_fun,
  cell_fun = cell_fun,
  cluster_rows = F,
  cluster_columns = T,
  row_split = c(rep("Skin from Male Mice", 4), rep("Skin from Female Mice", 6)),
  left_annotation = leftha,
  show_row_names = F,
  border = T,
  na_col = "grey",
  column_names_gp = gpar(fontsize = 12)
)

# Plot heatmap
pdf(file.path(output_dir, "pathway_enrichment_heatmap.pdf"))

draw(NES_heatmap, annotation_legend_list = significance_lgd)

dev.off()