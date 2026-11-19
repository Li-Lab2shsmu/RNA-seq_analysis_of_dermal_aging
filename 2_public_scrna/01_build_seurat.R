#!/usr/bin/env Rscript

# Purpose: build the Seurat object used as the BayesPrism reference.
# Input:   data/seurat_input/matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz,
#          metadata.csv, umap.csv, scvi.csv
# Output:  data/healthy_fibro_seurat.rds
# Usage:   Rscript 01_build_seurat.R
#
# Counts are read with Matrix::readMM and checked to be integers before the
# object is built.

suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
    stop("cannot locate this script, run it as: Rscript 01_build_seurat.R")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
setwd(script_dir)
source("config.R")

counts <- readMM(gzfile(file.path(SEURAT_DIR, "matrix.mtx.gz")))
genes <- read.table(gzfile(file.path(SEURAT_DIR, "features.tsv.gz")),
                    sep = "\t", stringsAsFactors = FALSE)[, 1]
cells <- read.table(gzfile(file.path(SEURAT_DIR, "barcodes.tsv.gz")),
                    sep = "\t", stringsAsFactors = FALSE)[, 1]

rownames(counts) <- genes
colnames(counts) <- cells
message(nrow(counts), " genes x ", ncol(counts), " cells")
message("integer counts: ", all(counts@x == floor(counts@x)))

meta <- read.csv(file.path(SEURAT_DIR, "metadata.csv"), row.names = 1,
                 stringsAsFactors = TRUE)

sobj <- CreateSeuratObject(counts = counts, assay = "RNA", meta.data = meta)
sobj <- NormalizeData(sobj, normalization.method = "LogNormalize",
                      scale.factor = 10000, verbose = FALSE)

umap_mat <- as.matrix(read.csv(file.path(SEURAT_DIR, "umap.csv"),
                               row.names = 1))
colnames(umap_mat) <- c("UMAP_1", "UMAP_2")
sobj[["umap"]] <- CreateDimReducObject(embeddings = umap_mat,
                                       key = "UMAP_", assay = "RNA")

scvi_df <- read.csv(file.path(SEURAT_DIR, "scvi.csv"), row.names = 1)
scvi_mat <- as.matrix(scvi_df)
colnames(scvi_mat) <- paste0("scVI_", seq_len(ncol(scvi_mat)))
sobj[["scvi"]] <- CreateDimReducObject(embeddings = scvi_mat,
                                       key = "scVI_", assay = "RNA")

saveRDS(sobj, file = SEURAT_RDS)
message("wrote ", SEURAT_RDS, " (",
        round(file.info(SEURAT_RDS)$size / 1e6, 1), " MB)")
