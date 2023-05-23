# Load libraries
library(reshape2)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(tidyr)
library(RColorBrewer)

# Set output directory
output_dir <- "plots"

########### PREPROCESSING ###########

# Load GSEA results from file
df_gsea <- read.csv("fgsea.csv", header = T)

# Create vector of all unique pathway names
pathway_names <- unique(df_gsea$pathway)

# create a list of all of the leading edge genes for each pathway
union_genes <- list()
for (p in pathway_names) {
  path_df <- df_gsea %>% dplyr::filter(pathway == p)
  genes <-
    unlist(strsplit(paste(path_df$leadingEdge, collapse = "|"), "\\|"))
  union_genes[[p]] <- unique(genes)
}

# Load deseq results file
deseq_df <- read.csv("deseq_annotated.csv", header = T)

# Reduce down to only HGNC mapped genes
deseq_df_mapped <- deseq_df[deseq_df$HGNC != "",]

# Return dataframes for just the leading edge genes for each pathway
pathway_dfs <-
  lapply(union_genes, function(x)
    deseq_df_mapped[deseq_df_mapped$HGNC %in% x,])

# Take a note of significant genes FDR <= 0.01
sig_gene_symbols <-
  lapply(pathway_dfs, function(x)
    subset(x, padj <= 0.01)$SYMBOL)

# Filter down to only genes that were significant in at least one subset
pathway_dfs_sig_genes <-
  mapply(
    function(x, y)
      subset(x, x$SYMBOL %in% y),
    x = pathway_dfs,
    y = sig_gene_symbols,
    SIMPLIFY = F
  )

# Ensure the subset order matches the original deseq_df
pathway_dfs_sig_genes <- lapply(pathway_dfs_sig_genes, function(x) {
  x %>% mutate(subset = factor(subset, levels = unique(deseq_df$subset)))
})

# Convert dataframes into t-score matrices
tscore_matrices <-
  lapply(pathway_dfs_sig_genes, function(x)
    t(acast(x, subset ~ SYMBOL, value.var = "stat")))

########### PLOTTING ###########

# Load metadata
meta_data_df <-
  read.csv("subsets.csv", row.names = 1, encoding = "UTF-8")

# Heatmap annotation bar colours
anno_colors <- list(
  Mission = c(
    "MHU-2" = "#1D91C0",
    "RR-5" = "#F7FCFD",
    "RR-7" = "#CB181D"
  ),
  Strain = c(
    "C57BL" = "#0072b2",
    "BALB" = "#914770",
    "C3H" = "#7f7f7f"
  ),
  Diet = c(
    "JC" = "#762A83",
    "JCwFOS" = "#FFEE99",
    "NuRFB" = "#225555"
  ),
  Euthanasia = c("LAR" = "#EE8866", "non-LAR" = "#EEDD88"),
  Duration = c("25" = "#f419ec", "30" = "#4419f4", "75" = "#19f4f2"),
  Tissue = c("Dorsal" = "#222255", "Femoral" = "#663333"),
  Age = c("9" = "#44BB99", "11" = "#99DDFF", "30" = "#77AADD"),
  Recovery = c("1" = "#F26E01", "30" = "#E34A27", "0" = "#FFFFFF")
)

# Heatmap annotation bar
leftha <- columnAnnotation(
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
  height = unit(4, "cm"),
  simple_anno_size_adjust = TRUE
)

# Legend for significance levels
significance_lgd <- Legend(
  labels = c("<= 0.05", "<= 0.01", "<= 0.001"),
  title = "Significance (FDR)",
  type = "points",
  pch = c("*", "**", "***"),
  background = "white"
)

# Function to add levels of significance
cell_fun <- function(j, i, x, y, width, height, fill) {
  padj_value <- deseq_df %>%
    dplyr::filter(SYMBOL == rownames(matrix)[i],
                  subset == colnames(matrix)[j]) %>%
    pull(padj)
  
  if (length(padj_value) == 0 || is.na(padj_value))
    return()
  
  y_adj <- y - unit(0.0075, "npc")
  
  if (padj_value <= 0.001)
    grid.text("***",
              x,
              y_adj,
              just = c("center", "bottom"),
              gp = gpar(fontsize = 12))
  else if (padj_value <= 0.01)
    grid.text("**",
              x,
              y_adj,
              just = c("center", "bottom"),
              gp = gpar(fontsize = 12))
  else if (padj_value <= 0.05)
    grid.text("*",
              x,
              y_adj,
              just = c("center", "bottom"),
              gp = gpar(fontsize = 12))
}

# Heatmap colour scale
col_fun <-
  colorRamp2(c(-15,-5, 0, 5, 15),
             c("darkblue", "blue", "white", "red", "darkred"))

# Plot heatmap into pdf file
pdf(file.path(output_dir, "gene_heatmaps.pdf"), height = 15)

text_size <-
  c(
    "Estrogen Signalling" = 10,
    "Insulin Resistance" = 10,
    "Insulin Signalling" = 7
  )
for (matrix_name in names(tscore_matrices)) {
  matrix <- tscore_matrices[[matrix_name]]
  
  tscore_heatmap <- Heatmap(
    matrix,
    name = "t-score",
    column_title = matrix_name,
    col = col_fun,
    cell_fun = cell_fun,
    cluster_rows = T,
    cluster_columns = F,
    column_split = c(rep("Male Mice", 4), rep("Female Mice", 6)),
    top_annotation = leftha,
    show_row_names = T,
    show_column_names = F,
    border = T,
    na_col = "grey",
    row_names_gp = gpar(fontsize = text_size[matrix_name])
  )
  
  draw(tscore_heatmap, annotation_legend_list = significance_lgd)
  
}
dev.off()