#!/usr/bin/env Rscript

# Purpose: compare TF regulatory influence between the Middle and Aged networks.
# Input:   results/genie3/GENIE3_hub_{Middle,Aged}.tsv
#          results/tf/TF_rank_table.tsv (optional)
# Output:  results/genie3/GENIE3_TF_outdegree.tsv
# Usage:   Rscript 14_genie3_tf_influence.R
#
# Influence is the TF out-degree after the top-k truncation. FCI is
# (Aged degree + eps) / (Middle degree + eps); a positive log2 FCI means the TF
# regulates more targets in Aged.

suppressPackageStartupMessages({
    library(dplyr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
    stop("cannot locate this script, run it as: Rscript 14_genie3_tf_influence.R")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
setwd(script_dir)
source("config.R")

input_dir <- file.path(DIR_RESULTS, "genie3")
tf_rank_file <- file.path(DIR_RESULTS, "tf", "TF_rank_table.tsv")
output_dir <- input_dir

eps <- 1

diff_tf <- character(0)
if (file.exists(tf_rank_file)) {
    tf_rank <- read.table(tf_rank_file, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)
    diff_tf <- tf_rank$GeneID[tf_rank[[SIG_CRITERIA]] < SIG_THRESHOLD]
}

read_hub <- function(net) {
    f <- file.path(input_dir, paste0("GENIE3_hub_", net, ".tsv"))
    if (!file.exists(f)) {
        stop("hub table not found: ", f,
             "\n  run 12_genie3_hub.R first")
    }
    h <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    d <- data.frame(TF = h$regulatoryGene, degree = h$out_degree,
                    stringsAsFactors = FALSE)
    colnames(d)[2] <- paste0("deg_", net)
    d
}

h_middle <- read_hub("Middle")
h_aged <- read_hub("Aged")

df <- merge(h_middle, h_aged, by = "TF", all = TRUE)
df[is.na(df)] <- 0

df$log2FCI <- log2((df$deg_Aged + eps) / (df$deg_Middle + eps))
df$abs_log2FCI <- abs(df$log2FCI)
df$is_diff_tf <- df$TF %in% diff_tf

out_table <- df %>% arrange(desc(abs_log2FCI))
write.table(out_table, file = file.path(output_dir, "GENIE3_TF_outdegree.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote ", file.path(output_dir, "GENIE3_TF_outdegree.tsv"),
        " (", nrow(out_table), " TFs)")
