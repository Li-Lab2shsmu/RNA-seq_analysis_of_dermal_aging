# Shared settings for the public scRNA-seq deconvolution pipeline.

# 10X-format export of the healthy skin fibroblast reference.
SEURAT_DIR <- "data/seurat_input"
SEURAT_RDS <- "data/healthy_fibro_seurat.rds"
RESULTS_DIR <- "results"
THETA_FILE <- file.path(RESULTS_DIR, "theta_cell_fractions.tsv")

# Bulk counts produced by 1_bulk_rnaseq/03_expression_matrix.R.
BULK_COUNT <- file.path("..", "1_bulk_rnaseq", "results", "MAT_COUNT_noD48.tsv")

# BayesPrism run settings.
N_CORES <- 8
