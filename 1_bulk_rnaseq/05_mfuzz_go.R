#!/usr/bin/env Rscript

# Purpose: GO enrichment for each Mfuzz cluster.
# Input:   results/mfuzz/Mfuzz_clustering_results_6.csv
# Output:  results/mfuzz_go/GO_analysis_clusters_6.xlsx
# Usage:   Rscript 05_mfuzz_go.R
#
# Only genes whose maximum membership exceeds the cutoff are tested, which
# keeps the core members of each cluster and drops genes shared between
# clusters.

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(dplyr)
    library(openxlsx)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
    stop("cannot locate this script, run it as: Rscript 05_mfuzz_go.R")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
setwd(script_dir)
source("config.R")

cluster_file <- file.path(DIR_RESULTS, "mfuzz", "Mfuzz_clustering_results_6.csv")
output_dir <- file.path(DIR_RESULTS, "mfuzz_go")

core_membership <- 0.45   # 0 keeps every gene
species_db <- "org.Hs.eg.db"
go_ont <- "ALL"
key_type <- "SYMBOL"
pvalue_cut <- 0.05
qvalue_cut <- 0.05

if (!file.exists(cluster_file)) {
    stop("cluster result not found: ", cluster_file,
         "\n  run 04_mfuzz_cluster.R first")
}
if (!is.numeric(core_membership) || core_membership < 0 || core_membership > 1) {
    stop("core_membership must be between 0 and 1, got: ", core_membership)
}

cluster_result <- read.csv(cluster_file, row.names = 1,
                           stringsAsFactors = FALSE)
n_clusters <- max(cluster_result$Cluster)
message(nrow(cluster_result), " genes in ", n_clusters, " clusters")

gene_cluster <- cluster_result$Cluster
names(gene_cluster) <- rownames(cluster_result)

if (core_membership > 0) {
    mem_cols <- paste0("X", seq_len(n_clusters))
    missing_cols <- setdiff(mem_cols, colnames(cluster_result))
    if (length(missing_cols) > 0) {
        stop("membership columns missing from the cluster result: ",
             paste(missing_cols, collapse = ", "))
    }

    membership_mat <- as.matrix(cluster_result[, mem_cols])
    max_membership <- apply(membership_mat, 1, max)
    core_flag <- max_membership > core_membership
    message("core genes: ", sum(core_flag), " / ", length(core_flag))

    gene_cluster <- gene_cluster[core_flag]
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

wb <- createWorkbook()
cluster_summary <- list()

for (k in seq_len(n_clusters)) {
    genes_k <- names(gene_cluster[gene_cluster == k])

    if (length(genes_k) < 5) {
        message("cluster ", k, ": ", length(genes_k),
                " genes, too few to test, skipped")
        next
    }

    ego <- enrichGO(
        gene          = genes_k,
        OrgDb         = species_db,
        keyType       = key_type,
        ont           = go_ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = pvalue_cut,
        qvalueCutoff  = qvalue_cut
    )

    if (is.null(ego) || nrow(ego) == 0) {
        message("cluster ", k, ": no enriched term")
        cluster_summary[[as.character(k)]] <- data.frame(
            Cluster = k, GeneCount = length(genes_k),
            GOTerms = 0, stringsAsFactors = FALSE
        )
        next
    }

    ego_result <- ego@result %>% arrange(qvalue)
    bp_n <- sum(ego_result$ONTOLOGY == "BP")
    cc_n <- sum(ego_result$ONTOLOGY == "CC")
    mf_n <- sum(ego_result$ONTOLOGY == "MF")

    sheet_name <- paste0("Cluster", k)
    addWorksheet(wb, sheetName = sheet_name)
    writeData(wb, sheet = sheet_name, x = ego_result)
    setColWidths(wb, sheet_name, cols = seq_len(ncol(ego_result)), widths = "auto")

    cluster_summary[[as.character(k)]] <- data.frame(
        Cluster = k,
        GeneCount = length(genes_k),
        BP = bp_n, CC = cc_n, MF = mf_n, Total = nrow(ego_result),
        stringsAsFactors = FALSE
    )

    message("cluster ", k, ": ", nrow(ego_result), " terms")
}

if (length(cluster_summary) > 0) {
    summary_df <- do.call(rbind, cluster_summary)
    addWorksheet(wb, "Summary", gridLines = TRUE)
    writeData(wb, "Summary", x = summary_df)
    setColWidths(wb, "Summary", cols = seq_len(ncol(summary_df)), widths = "auto")
}

xlsx_file <- file.path(output_dir, "GO_analysis_clusters_6.xlsx")
saveWorkbook(wb, xlsx_file, overwrite = TRUE)
message("wrote ", xlsx_file)
