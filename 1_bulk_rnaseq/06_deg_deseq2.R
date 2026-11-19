#!/usr/bin/env Rscript

# Purpose: DESeq2 differential expression of Aged against Middle.
# Input:   results/MAT_COUNT_noD48.tsv
# Output:  results/deg/Aged_vs_Middle_results.csv, results/deg/all_results.xlsx
# Usage:   Rscript 06_deg_deseq2.R
#
# Low-count genes are removed first. log2FoldChange > 0 means higher in Aged.

suppressPackageStartupMessages({
    library(DESeq2)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
    stop("cannot locate this script, run it as: Rscript 06_deg_deseq2.R")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
setwd(script_dir)
source("config.R")

counts_file <- MAT_COUNT
output_dir <- file.path(DIR_RESULTS, "deg")
comparison_name <- paste0(DEG_TEST, "_vs_", DEG_CONTROL)

if (!SIG_CRITERIA %in% c("padj", "pvalue")) {
    stop("SIG_CRITERIA must be padj or pvalue, got: ", SIG_CRITERIA)
}
if (!file.exists(counts_file)) {
    stop("count matrix not found: ", counts_file,
         "\n  run 03_expression_matrix.R first")
}

counts <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t",
                     stringsAsFactors = FALSE, check.names = FALSE)

# Keep genes with at least MIN_COUNT reads in one sample.
keep <- rowSums(counts >= MIN_COUNT) >= 1
counts <- counts[keep, , drop = FALSE]
counts <- counts[rowSums(counts) != 0, , drop = FALSE]
if (nrow(counts) == 0) {
    stop("no gene left after low-count filtering")
}

test_samples <- GROUPS[[DEG_TEST]]
control_samples <- GROUPS[[DEG_CONTROL]]
selected <- c(test_samples, control_samples)
missing_samples <- setdiff(selected, colnames(counts))
if (length(missing_samples) > 0) {
    stop("samples missing from the count matrix: ",
         paste(missing_samples, collapse = ", "))
}

counts <- counts[, selected, drop = FALSE]

condition <- factor(c(rep(DEG_TEST, length(test_samples)),
                      rep(DEG_CONTROL, length(control_samples))),
                    levels = c(DEG_TEST, DEG_CONTROL))
col_data <- data.frame(row.names = colnames(counts), condition = condition)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = col_data,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", DEG_TEST, DEG_CONTROL))

res_df <- as.data.frame(res)
res_df$GeneID <- rownames(res_df)
res_df <- res_df[, c("GeneID", "baseMean", "log2FoldChange", "lfcSE",
                     "stat", "pvalue", "padj")]
res_df <- res_df[order(res_df[, SIG_CRITERIA], na.last = TRUE), ]
res_df$Significant <- ifelse(
    res_df[, SIG_CRITERIA] < SIG_THRESHOLD & !is.na(res_df[, SIG_CRITERIA]),
    "Yes", "No"
)

n_sig <- sum(res_df[, SIG_CRITERIA] < SIG_THRESHOLD, na.rm = TRUE)
n_up <- sum(res_df[, SIG_CRITERIA] < SIG_THRESHOLD &
            res_df$log2FoldChange > 0, na.rm = TRUE)
n_down <- sum(res_df[, SIG_CRITERIA] < SIG_THRESHOLD &
              res_df$log2FoldChange < 0, na.rm = TRUE)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
csv_file <- file.path(output_dir, paste0(comparison_name, "_results.csv"))
write.csv(res_df, file = csv_file, row.names = FALSE)
message("wrote ", csv_file)
message(n_sig, " genes with ", SIG_CRITERIA, " < ", SIG_THRESHOLD,
        " (up ", n_up, ", down ", n_down, ")")

if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, comparison_name)
    openxlsx::writeData(wb, comparison_name, res_df)
    openxlsx::setColWidths(wb, comparison_name,
                           cols = seq_len(ncol(res_df)), widths = "auto")
    xlsx_file <- file.path(output_dir, "all_results.xlsx")
    openxlsx::saveWorkbook(wb, xlsx_file, overwrite = TRUE)
    message("wrote ", xlsx_file)
}
