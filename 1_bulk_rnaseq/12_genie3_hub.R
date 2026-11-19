#!/usr/bin/env Rscript

# Purpose: sparsify each GENIE3 network and identify hub transcription factors.
# Input:   results/genie3/GENIE3_linklist_{Middle,Aged}.tsv
#          results/tf/TF_rank_table.tsv (optional, for the differential flag)
# Output:  results/genie3/GENIE3_hub_{Middle,Aged}.tsv
#          results/genie3/GENIE3_hub_all.tsv
# Usage:   Rscript 12_genie3_hub.R
#
# GENIE3 weights are only a ranking, so each target keeps its top_k regulators
# before a TF's out-degree is counted.

suppressPackageStartupMessages({
    library(dplyr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
    stop("cannot locate this script, run it as: Rscript 12_genie3_hub.R")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
setwd(script_dir)
source("config.R")

input_dir <- file.path(DIR_RESULTS, "genie3")
tf_rank_file <- file.path(DIR_RESULTS, "tf", "TF_rank_table.tsv")
output_dir <- input_dir

top_k <- 10

diff_tf <- NULL
if (file.exists(tf_rank_file)) {
    tf_rank <- read.table(tf_rank_file, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)
    diff_tf <- tf_rank$GeneID[tf_rank[[SIG_CRITERIA]] < SIG_THRESHOLD]
    message(length(diff_tf), " differential TFs")
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
hub_summary <- list()

for (net in GENIE3_NETWORKS) {
    link_file <- file.path(input_dir, paste0("GENIE3_linklist_", net, ".tsv"))
    if (!file.exists(link_file)) {
        stop("link list not found: ", link_file,
             "\n  run 11_genie3_network.R first")
    }

    link <- read.table(link_file, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)

    sparse <- link %>%
        group_by(targetGene) %>%
        slice_max(weight, n = top_k, with_ties = TRUE) %>%
        ungroup()

    hub <- sparse %>%
        group_by(regulatoryGene) %>%
        summarise(
            out_degree  = n_distinct(targetGene),
            mean_weight = mean(weight),
            .groups     = "drop"
        ) %>%
        arrange(desc(out_degree), desc(mean_weight))
    hub$is_diff_tf <- hub$regulatoryGene %in% diff_tf

    hub_file <- file.path(output_dir, paste0("GENIE3_hub_", net, ".tsv"))
    write.table(hub, file = hub_file, sep = "\t", quote = FALSE, row.names = FALSE)
    message(net, ": ", nrow(sparse), " edges, ", nrow(hub), " TFs")

    hub_summary[[net]] <- hub %>% mutate(network = net)
}

all_hub <- bind_rows(hub_summary)
summary_file <- file.path(output_dir, "GENIE3_hub_all.tsv")
write.table(all_hub, file = summary_file, sep = "\t", quote = FALSE,
            row.names = FALSE)
message("wrote ", summary_file)
