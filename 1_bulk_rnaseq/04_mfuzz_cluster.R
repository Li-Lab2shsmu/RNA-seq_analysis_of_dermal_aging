#!/usr/bin/env Rscript

# Purpose: Mfuzz soft clustering of age-group mean expression profiles.
# Input:   results/MAT_RPKM_noD48.tsv
# Output:  results/mfuzz/Mfuzz_clustering_results_<k>.csv
# Usage:   Rscript 04_mfuzz_cluster.R
#
# Clustering runs on the mean profile of each age group, so the four group
# means are the time axis. Membership and cluster assignment for every gene are
# written for each candidate cluster count.

suppressPackageStartupMessages({
    library(dplyr)
    library(Mfuzz)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
    stop("cannot locate this script, run it as: Rscript 04_mfuzz_cluster.R")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
setwd(script_dir)
source("config.R")

input_file <- MAT_RPKM
output_dir <- file.path(DIR_RESULTS, "mfuzz")

# Candidate cluster counts. k = 6 is the one used in the manuscript.
cluster_numbers <- c(5, 6, 7, 8, 9, 10)
seed_value <- 123
min_std <- 0.1

if (!file.exists(input_file)) {
    stop("input matrix not found: ", input_file,
         "\n  run 03_expression_matrix.R first")
}

message("reading ", input_file)
data <- read.table(input_file, header = TRUE, row.names = 1,
                   sep = "\t", check.names = FALSE)
message(nrow(data), " genes x ", ncol(data), " samples")

group_cols <- unlist(GROUPS, use.names = FALSE)
missing <- setdiff(group_cols, colnames(data))
if (length(missing) > 0) {
    stop("group samples not found in the matrix: ",
         paste(missing, collapse = ", "))
}

# The group means are the time axis, in the order defined by GROUPS.
group_means <- lapply(GROUPS, function(samples) {
    rowMeans(data[, samples, drop = FALSE], na.rm = TRUE)
})
expr_matrix <- as.matrix(as.data.frame(group_means))
colnames(expr_matrix) <- names(GROUPS)

message("group-mean matrix: ", nrow(expr_matrix), " genes x ",
        ncol(expr_matrix), " time points")

eset <- new("ExpressionSet", exprs = expr_matrix)
eset <- fill.NA(eset, mode = "mean")
eset <- filter.std(eset, min.std = min_std, visu = FALSE)
eset <- standardise(eset)
message("genes after the standard-deviation filter: ", nrow(exprs(eset)))

m_value <- mestimate(eset)
message("estimated fuzzifier m: ", round(m_value, 6))

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
time_labels <- names(GROUPS)

for (k in cluster_numbers) {
    message("cluster count C = ", k)

    set.seed(seed_value)
    clusters <- mfuzz(eset, c = k, m = m_value)

    cluster_result <- data.frame(
        GeneID = rownames(clusters$membership),
        Cluster = clusters$cluster,
        clusters$membership,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    rownames(cluster_result) <- NULL

    csv_file <- file.path(output_dir,
                          paste0("Mfuzz_clustering_results_", k, ".csv"))
    write.csv(cluster_result, file = csv_file, row.names = FALSE)
    message("  wrote ", csv_file)
}

message("done")
