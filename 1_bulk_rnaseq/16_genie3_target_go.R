#!/usr/bin/env Rscript

# Purpose: GO enrichment of the top GENIE3 targets of NR4A1, MYC and ETV4 in
#          the Middle network.
# Input:   results/genie3/GENIE3_linklist_Middle.tsv
# Output:  results/genie3_go/<TF>_target_GO.csv
#          results/genie3_go/all_target_GO.csv
# Usage:   Rscript 16_genie3_target_go.R
#
# The three TFs are down-regulated with age, so their Middle-network targets
# describe the regulatory programme that is lost.

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(dplyr)
})

options(warn = -1)
on.exit(options(warn = 0))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
    stop("cannot locate this script, run it as: Rscript 16_genie3_target_go.R")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
setwd(script_dir)
source("config.R")

link_file <- file.path(DIR_RESULTS, "genie3", "GENIE3_linklist_Middle.tsv")
output_dir <- file.path(DIR_RESULTS, "genie3_go")

tfs <- c("NR4A1", "MYC", "ETV4")
top_targets <- 100

species_db <- "org.Hs.eg.db"
key_type <- "SYMBOL"
ont <- "BP"
pvalue_cut <- 0.05
qvalue_cut <- 0.2

if (!file.exists(link_file)) {
    stop("link list not found: ", link_file,
         "\n  run 11_genie3_network.R first")
}

link <- read.table(link_file, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

target_sets <- list()
for (tf in tfs) {
    sub <- link[link$regulatoryGene == tf, , drop = FALSE]
    sub <- sub[order(-sub$weight), , drop = FALSE]
    target_sets[[tf]] <- head(unique(sub$targetGene), top_targets)
    message(tf, ": ", length(target_sets[[tf]]), " targets")
}

all_results <- list()

for (tf in tfs) {
    genes <- target_sets[[tf]]
    if (length(genes) < 5) {
        message(tf, ": fewer than 5 targets, skipped")
        next
    }

    ego <- enrichGO(
        gene          = genes,
        OrgDb         = species_db,
        keyType       = key_type,
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = pvalue_cut,
        qvalueCutoff  = qvalue_cut
    )
    if (is.null(ego) || nrow(ego) == 0) {
        message(tf, ": no enriched term")
        next
    }

    res_df <- as.data.frame(ego) %>% arrange(p.adjust)
    res_df$TF <- tf
    all_results[[tf]] <- res_df

    write.csv(res_df,
              file = file.path(output_dir, paste0(tf, "_target_GO.csv")),
              row.names = FALSE)
}

if (length(all_results) > 0) {
    combined <- do.call(rbind, all_results)
    write.csv(combined, file = file.path(output_dir, "all_target_GO.csv"),
              row.names = FALSE)
    message("wrote ", file.path(output_dir, "all_target_GO.csv"))
}
