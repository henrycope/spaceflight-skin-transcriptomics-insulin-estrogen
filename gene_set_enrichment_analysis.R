# Libraries & file/dir names
library("fgsea")
library("dplyr")
library("data.table")

deseq_in_fn <- "deseq.csv"
gmt_dir <- "gmt_files"
gsea_out_fn <- "fgsea.csv"

set.seed(42)

# Check input data exists
if (!file.exists(deseq_in_fn) || !dir.exists(gmt_dir)) {
  stop("Input file or directory not found")
}

# Load gmt files
gene_sets <-
  lapply(list.files(gmt_dir, full.names = T),
         fgsea::gmtPathways)
names(gene_sets) <- list.files(gmt_dir)

# Load deseq results file
deseq_df <- read.csv(deseq_in_fn, header = T)

# Drop all columns other than HGNC, t-score and subset
deseq_df <- deseq_df[, c("HGNC", "stat", "subset")]

# Remove unmapped genes
deseq_df <- deseq_df[deseq_df$HGNC != "",]

# Split into seperate dfs on subset
deseq_dfs <-
  split(deseq_df,
        factor(deseq_df$subset, levels = (unique(
          deseq_df$subset
        ))))

# Sort by t-score, descending
sorted_deseq_dfs <-
  lapply(deseq_dfs, function(x)
    x[order(-x$stat), c("HGNC", "stat")])

# Remove duplicates
unique_deseq_dfs <-
  lapply(sorted_deseq_dfs, function(x)
    x[!duplicated(x$HGNC), ])

# Create t-score rank vectors
rank_vectors <-
  lapply(unique_deseq_dfs, function(x)
    setNames(x$stat, x$HGNC))

# Perform GSEA on rank vectors (all combos of rank vectors and pathways)
# TODO: remove nested lapply, expand.grid?
fgsea_res <- lapply(gene_sets,
                    function(x)
                      lapply(rank_vectors, function(y)
                        fgsea(
                          pathways = x,
                          stats = y,
                          minSize = 1,
                          maxSize = 10000
                        )))

# Condense to a combined df per category
category_dfs <- lapply(fgsea_res, function(x) dplyr::bind_rows(x, .id = "subset"))

# Condense to one combined df
combined_df <- dplyr::bind_rows(category_dfs, .id = "category")
combined_df$category <- gsub( ".gmt", "", combined_df$category)

# Save combined df
fwrite(combined_df, gsea_out_fn)