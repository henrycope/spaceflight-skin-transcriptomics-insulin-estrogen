library(biomaRt)
library(dplyr)
library(tidyr)

# Read the csv file
data <- read.csv("deseq.csv")

get_gene_symbols <- function(ensembl_ids, mart) {
  # Retrieve gene symbols
  gene_symbols <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = mart
  )
  return(gene_symbols)
}

mouse_ensembl_ids <- unique(data$ENSEMBL)
mouse_ensembl_ids <- mouse_ensembl_ids[!is.na(mouse_ensembl_ids)]

# Get mouse gene symbols
mart_mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
mouse_gene_symbols <- get_gene_symbols(mouse_ensembl_ids, mart_mouse)

# Get human gene symbols by first converting the mouse IDs to human
mart_human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
human_orthologs <- getBM(
  attributes = c('ensembl_gene_id', 'hsapiens_homolog_ensembl_gene', 'hsapiens_homolog_perc_id'),
  filters = 'ensembl_gene_id',
  values = mouse_ensembl_ids,
  mart = mart_mouse
)

# Select the human ortholog with the highest homology score for each mouse gene
human_orthologs <- human_orthologs %>% group_by(ensembl_gene_id) %>% top_n(1, hsapiens_homolog_perc_id)

# Get HGNC symbol for human ENSEMBL ids - first column is human ensembl id
human_gene_symbols <- get_gene_symbols(human_orthologs$hsapiens_homolog_ensembl_gene, mart_human)

# Connect human symbols to mice ensembl ids, and remove duplicates
human_gene_symbols <- left_join(human_orthologs, human_gene_symbols, by = c("hsapiens_homolog_ensembl_gene" = "ensembl_gene_id"))
human_gene_symbols <- human_gene_symbols %>% distinct(ensembl_gene_id, .keep_all = TRUE)

# Merge mouse gene symbols with original data
data_with_mouse_symbols <- left_join(data, mouse_gene_symbols, by = c("ENSEMBL" = "ensembl_gene_id"))

# Merge human gene symbols with the data
data_with_both_symbols <- left_join(data_with_mouse_symbols, human_gene_symbols, c("ENSEMBL" = "ensembl_gene_id"))

# Rename the columns
colnames(data_with_both_symbols)[colnames(data_with_both_symbols) == "external_gene_name.x"] <- "SYMBOL"
colnames(data_with_both_symbols)[colnames(data_with_both_symbols) == "external_gene_name.y"] <- "HGNC"

# Replace NAs with blanks
data_with_both_symbols[] <- lapply(data_with_both_symbols, function(x) ifelse(is.na(x), "", x))

# Write the output to a new csv file
write.csv(data_with_both_symbols, "deseq_annotated.csv", row.names = FALSE)