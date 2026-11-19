#!/usr/bin/env Rscript

# Purpose: compare the regulatory position of TFs between the Middle and Aged
#          networks and write the differential hubs.
# Input:   results/genie3/GENIE3_hub_{Middle,Aged}.tsv
#          results/tf/TF_rank_table.tsv (optional)
# Output:  results/genie3/GENIE3_diffhub_Middle_Aged.tsv
# Usage:   Rscript 15_genie3_diffhub.R
#
# The two networks have different out-degree ranges, so position is compared as
# a within-network percentile rank. delta is percentile_Aged - percentile_Middle.

suppressPackageStartupMessages({
    library(dplyr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
    stop("cannot locate this script, run it as: Rscript 15_genie3_diffhub.R")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
setwd(script_dir)
source("config.R")

input_dir <- file.path(DIR_RESULTS, "genie3")
tf_rank_file <- file.path(DIR_RESULTS, "tf", "TF_rank_table.tsv")
output_dir <- input_dir

# TFs below this out-degree in either network are too unstable to compare.
min_degree <- 20

read_hub <- function(net) {
    f <- file.path(input_dir, paste0("GENIE3_hub_", net, ".tsv"))
    if (!file.exists(f)) {
        stop("hub table not found: ", f,
             "\n  run 12_genie3_hub.R first")
    }
    read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}

middle <- read_hub("Middle")
aged <- read_hub("Aged")

middle$percentile <- rank(middle$out_degree, ties.method = "average") / nrow(middle)
aged$percentile <- rank(aged$out_degree, ties.method = "average") / nrow(aged)

df <- merge(
    middle[, c("regulatoryGene", "out_degree", "percentile")],
    aged[, c("regulatoryGene", "out_degree", "percentile")],
    by = "regulatoryGene", all = TRUE, suffixes = c("_Middle", "_Aged")
)
df[is.na(df)] <- 0

df$delta <- df$percentile_Aged - df$percentile_Middle

diff_tf <- character(0)
if (file.exists(tf_rank_file)) {
    tf_rank <- read.table(tf_rank_file, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)
    diff_tf <- tf_rank$GeneID[tf_rank[[SIG_CRITERIA]] < SIG_THRESHOLD]
}
df$is_diff_tf <- df$regulatoryGene %in% diff_tf

df <- df[df$out_degree_Middle >= min_degree &
         df$out_degree_Aged >= min_degree, , drop = FALSE]
df <- df[order(-abs(df$delta)), , drop = FALSE]

out_table <- df[, c("regulatoryGene", "out_degree_Middle", "out_degree_Aged",
                    "percentile_Middle", "percentile_Aged", "delta",
                    "is_diff_tf")]
colnames(out_table) <- c("TF", "out_degree_Middle", "out_degree_Aged",
                         "percentile_Middle", "percentile_Aged", "delta",
                         "is_diff_tf")

write.table(out_table,
            file = file.path(output_dir, "GENIE3_diffhub_Middle_Aged.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote ", file.path(output_dir, "GENIE3_diffhub_Middle_Aged.tsv"),
        " (", nrow(out_table), " TFs)")
