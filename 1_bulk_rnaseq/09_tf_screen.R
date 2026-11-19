#!/usr/bin/env Rscript

# Purpose: extract transcription factors from the Aged vs Middle DEG table.
# Input:   results/deg/Aged_vs_Middle_results.csv
# Output:  results/tf/TF_rank_table.tsv
# Usage:   Rscript 09_tf_screen.R
#
# A gene counts as a transcription factor when it is annotated with GO:0003700
# or GO:0000981. Only TFs below the significance cutoff are kept; duplicates
# keep the lowest p-value.

suppressPackageStartupMessages({
    library(org.Hs.eg.db)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
    stop("cannot locate this script, run it as: Rscript 09_tf_screen.R")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
setwd(script_dir)
source("config.R")

deg_file <- file.path(DIR_RESULTS, "deg",
                      paste0(DEG_TEST, "_vs_", DEG_CONTROL, "_results.csv"))
output_dir <- file.path(DIR_RESULTS, "tf")

if (!file.exists(deg_file)) {
    stop("DEG result not found: ", deg_file,
         "\n  run 06_deg_deseq2.R first")
}

tf_ids <- unique(c(
    AnnotationDbi::select(org.Hs.eg.db,
                          keys = "GO:0003700",
                          columns = "SYMBOL",
                          keytype = "GOALL")$SYMBOL,
    AnnotationDbi::select(org.Hs.eg.db,
                          keys = "GO:0000981",
                          columns = "SYMBOL",
                          keytype = "GOALL")$SYMBOL
))
tf_set <- unique(tf_ids[!is.na(tf_ids)])
message(length(tf_set), " annotated transcription factors")

deg <- read.csv(deg_file, header = TRUE, stringsAsFactors = FALSE)
deg_tf <- deg[deg$GeneID %in% tf_set, , drop = FALSE]
deg_tf <- deg_tf[!is.na(deg_tf[[SIG_CRITERIA]]) &
                 deg_tf[[SIG_CRITERIA]] < SIG_THRESHOLD, , drop = FALSE]

deg_tf <- deg_tf[order(deg_tf$pvalue), ]
deg_tf <- deg_tf[!duplicated(deg_tf$GeneID), ]

up <- sum(deg_tf$log2FoldChange > 0)
down <- sum(deg_tf$log2FoldChange < 0)
message(nrow(deg_tf), " TFs pass the cutoff (up ", up, ", down ", down, ")")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
tf_table <- deg_tf[, c("GeneID", "log2FoldChange", "pvalue", "padj", "baseMean")]
tf_table <- tf_table[order(tf_table$pvalue), ]
write.table(tf_table,
            file = file.path(output_dir, "TF_rank_table.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote ", file.path(output_dir, "TF_rank_table.tsv"))
