#!/usr/bin/env Rscript

# Purpose: deconvolve the bulk samples with BayesPrism against the healthy
#          fibroblast reference and export the cell-type proportions.
# Input:   data/healthy_fibro_seurat.rds
#          ../1_bulk_rnaseq/results/MAT_COUNT_noD48.tsv
# Output:  results/theta_cell_fractions.tsv
# Usage:   Rscript 02_run_bayesprism.R
#
# BayesPrism needs a cell x gene reference and a sample x gene mixture, so both
# matrices are transposed after intersecting genes. The Gibbs sampler is
# seeded indirectly by the fixed chain settings; the run is long, keep N_CORES.

suppressPackageStartupMessages({
    library(BayesPrism)
    library(Seurat)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
    stop("cannot locate this script, run it as: Rscript 02_run_bayesprism.R")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
setwd(script_dir)
source("config.R")

if (!file.exists(SEURAT_RDS)) {
    stop("Seurat reference not found: ", SEURAT_RDS,
         "\n  run 01_build_seurat.R first")
}
if (!file.exists(BULK_COUNT)) {
    stop("bulk count matrix not found: ", BULK_COUNT,
         "\n  run 1_bulk_rnaseq/03_expression_matrix.R first")
}

sobj <- readRDS(SEURAT_RDS)
sc_counts <- GetAssayData(sobj, slot = "counts", assay = "RNA")
cell_labels <- as.character(sobj$celltype)
message(ncol(sobj), " cells x ", nrow(sobj), " genes, ",
        length(unique(cell_labels)), " cell types")

bulk <- read.table(BULK_COUNT, header = TRUE, row.names = 1, sep = "\t",
                   check.names = FALSE, stringsAsFactors = FALSE)
bulk <- as.matrix(bulk)
message(nrow(bulk), " bulk genes x ", ncol(bulk), " samples: ",
        paste(colnames(bulk), collapse = " "))

common_genes <- intersect(rownames(sc_counts), rownames(bulk))
if (length(common_genes) < 1000) {
    stop("only ", length(common_genes),
         " genes in common, check that both matrices use gene symbols")
}
message(length(common_genes), " genes in common")

sc_filt_t <- t(sc_counts[common_genes, ])
bulk_filt_t <- t(bulk[common_genes, , drop = FALSE])
rm(sc_counts, bulk)
gc()

n_cores <- max(1L, min(N_CORES, parallel::detectCores() - 1L))
message("running BayesPrism on ", n_cores, " cores")

my_prism <- new.prism(
    reference         = sc_filt_t,
    mixture           = bulk_filt_t,
    key               = NULL,
    cell.type.labels  = cell_labels,
    cell.state.labels = cell_labels,
    outlier.cut       = 0.01,
    outlier.fraction  = 0.1
)
rm(sc_filt_t, bulk_filt_t)
gc()

bp_res <- run.prism(
    prism    = my_prism,
    n.cores  = n_cores,
    update.gibbs = TRUE,
    gibbs.control = list(
        chain.length = 1000,
        burn.in      = 500,
        thinning     = 1
    )
)

theta <- bp_res@posterior.theta_f@theta
message(nrow(theta), " cell types x ", ncol(theta), " samples")

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
write.table(round(theta, 4), file = THETA_FILE, sep = "\t", quote = FALSE,
            col.names = NA)
message("wrote ", THETA_FILE)
