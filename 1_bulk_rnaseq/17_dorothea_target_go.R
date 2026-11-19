#!/usr/bin/env Rscript

# Purpose: GO enrichment of DoRothEA high-confidence targets of NR4A1, MYC and
#          ETV4.
# Input:   DoRothEA regulons (levels A and B)
# Output:  results/genie3_go/<TF>_dorothea_GO.csv
#          results/genie3_go/all_dorothea_GO.csv
# Usage:   Rscript 17_dorothea_target_go.R
#
# DoRothEA regulons are literature-curated, so they complement the inferred
# GENIE3 targets with higher-confidence TF-target links.

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(dplyr)
})

options(warn = -1)
on.exit(options(warn = 0))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
    stop("cannot locate this script, run it as: Rscript 17_dorothea_target_go.R")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
setwd(script_dir)
source("config.R")

output_dir <- file.path(DIR_RESULTS, "genie3_go")
tfs <- c("NR4A1", "MYC", "ETV4")

species_db <- "org.Hs.eg.db"
key_type <- "SYMBOL"
ont <- "BP"
pvalue_cut <- 0.05
qvalue_cut <- 0.2

regulons <- tryCatch(
    decoupleR::get_dorothea(organism = "human", levels = c("A", "B")),
    error = function(e) {
        message("online CollecTRI fetch failed, using bundled dorothea_hs")
        data("dorothea_hs", package = "dorothea")
        d <- dorothea_hs
        d[d$confidence %in% c("A", "B"), , drop = FALSE]
    }
)

if (!all(c("tf", "target") %in% colnames(regulons))) {
    stop("cannot identify tf/target columns in the regulon table")
}
regulons <- regulons[, c("tf", "target"), drop = FALSE]

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
all_results <- list()

for (tf in tfs) {
    genes <- unique(regulons$target[regulons$tf == tf])
    message(tf, ": ", length(genes), " DoRothEA targets")

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
              file = file.path(output_dir, paste0(tf, "_dorothea_GO.csv")),
              row.names = FALSE)
}

if (length(all_results) > 0) {
    combined <- do.call(rbind, all_results)
    write.csv(combined, file = file.path(output_dir, "all_dorothea_GO.csv"),
              row.names = FALSE)
    message("wrote ", file.path(output_dir, "all_dorothea_GO.csv"))
}
