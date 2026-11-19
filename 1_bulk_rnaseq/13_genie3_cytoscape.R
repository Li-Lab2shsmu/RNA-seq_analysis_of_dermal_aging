#!/usr/bin/env Rscript

# Purpose: export each GENIE3 network as Cytoscape edge and node tables.
# Input:   results/genie3/GENIE3_linklist_{Middle,Aged}.tsv
#          results/tf/TF_rank_table.tsv (optional)
# Output:  results/genie3/GENIE3_cytoscape_<network>_edges.tsv
#          results/genie3/GENIE3_cytoscape_<network>_nodes.tsv
# Usage:   Rscript 13_genie3_cytoscape.R
#
# Each target keeps its three strongest regulators so the exported network stays
# connected enough to be readable.

suppressPackageStartupMessages({
    library(dplyr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
    stop("cannot locate this script, run it as: Rscript 13_genie3_cytoscape.R")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
setwd(script_dir)
source("config.R")

input_dir <- file.path(DIR_RESULTS, "genie3")
tf_rank_file <- file.path(DIR_RESULTS, "tf", "TF_rank_table.tsv")
output_dir <- input_dir

# TFs discussed in the manuscript, flagged in the node table.
core_tfs <- c("BCL11A", "NR1D1", "RXRG", "ZBTB18",
              "TFAP2A", "JUN", "FOXQ1", "RFX8")

diff_tf <- character(0)
if (file.exists(tf_rank_file)) {
    tf_rank <- read.table(tf_rank_file, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)
    diff_tf <- tf_rank$GeneID[tf_rank[[SIG_CRITERIA]] < SIG_THRESHOLD]
}

for (net in GENIE3_NETWORKS) {
    link_file <- file.path(input_dir, paste0("GENIE3_linklist_", net, ".tsv"))
    if (!file.exists(link_file)) {
        stop("link list not found: ", link_file,
             "\n  run 11_genie3_network.R first")
    }

    link <- read.table(link_file, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)

    edges <- link %>%
        group_by(targetGene) %>%
        arrange(desc(weight)) %>%
        slice(1:3) %>%
        ungroup() %>%
        dplyr::select(source = regulatoryGene, target = targetGene, weight)

    tfs <- unique(edges$source)
    targets <- unique(edges$target)
    all_nodes <- unique(c(tfs, targets))
    out_deg <- edges %>% count(source, name = "out_degree")

    node_df <- data.frame(
        gene = all_nodes,
        node_type = ifelse(all_nodes %in% tfs, "TF", "target"),
        stringsAsFactors = FALSE
    ) %>%
        left_join(out_deg, by = c("gene" = "source")) %>%
        mutate(
            out_degree = ifelse(is.na(out_degree), 0, out_degree),
            is_diff_tf = gene %in% diff_tf,
            is_core = gene %in% core_tfs
        )

    edge_file <- file.path(output_dir,
                           paste0("GENIE3_cytoscape_", net, "_edges.tsv"))
    node_file <- file.path(output_dir,
                           paste0("GENIE3_cytoscape_", net, "_nodes.tsv"))
    write.table(edges, file = edge_file, sep = "\t", quote = FALSE,
                row.names = FALSE)
    write.table(node_df, file = node_file, sep = "\t", quote = FALSE,
                row.names = FALSE)
    message(net, ": ", nrow(edges), " edges, ", nrow(node_df), " nodes")
}
