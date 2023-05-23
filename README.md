# spaceflight-skin-transcriptomics-insulin-estrogen

Execution Order:
1) download_datasets.R - downloads datasets from GeneLab
2) deseq_pipeline.Rmd - Performs differential gene expression analysis from raw counts data
3) symbol_mapping.R - Maps ENSEMBL gene ids to mouse and human gene symbols
4) gene_set_enrichment_analysis.R - perform gene set enrichment analysis
5) pathway_enrichment_heatmap.R - produce a heatmap showing enrichment values for each pathway within each data subset
6) gene_level_plot.R - produce a set of heatmaps (one per pathway), showing the most significant genes from the pathwway
