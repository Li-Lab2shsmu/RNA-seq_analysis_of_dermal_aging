#!/usr/bin/env Rscript

# Purpose: GO enrichment of up- and down-regulated genes from the Aged vs
#          Middle comparison, run once on pvalue and once on padj.
# Input:   results/deg/Aged_vs_Middle_results.csv
# Output:  results/deg_go/Aged_vs_Middle_GO_enrichment_pvalue.csv and .xlsx
#          results/deg_go/Aged_vs_Middle_GO_enrichment_padj.csv and .xlsx
# Usage:   Rscript 07_deg_go.R
#
# Up- and down-regulated genes are enriched separately, and each round of
# filtering writes its own pair of files.

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(openxlsx)
})

options(warn = -1)
on.exit(options(warn = 0))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
    stop("cannot locate this script, run it as: Rscript 07_deg_go.R")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
setwd(script_dir)
source("config.R")

input_file <- file.path(DIR_RESULTS, "deg",
                        paste0(DEG_TEST, "_vs_", DEG_CONTROL, "_results.csv"))
output_dir <- file.path(DIR_RESULTS, "deg_go")
comp_name <- paste0(DEG_TEST, "_vs_", DEG_CONTROL)

sig_rounds <- c("pvalue", "padj")
sig_threshold <- 0.05

species_db <- "org.Hs.eg.db"
go_ont <- "ALL"
key_type <- "SYMBOL"
pvalue_cut <- 0.05
qvalue_cut <- 0.2

if (!file.exists(input_file)) {
    stop("DEG result not found: ", input_file,
         "\n  run 06_deg_deseq2.R first")
}

ont_list <- if (go_ont == "ALL") c("BP", "MF", "CC") else go_ont

run_enrich_go <- function(gene_list, ont) {
    if (length(gene_list) < 5) {
        message("    fewer than 5 genes, skipped")
        return(NULL)
    }
    ego <- enrichGO(
        gene          = gene_list,
        OrgDb         = species_db,
        keyType       = key_type,
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = pvalue_cut,
        qvalueCutoff  = qvalue_cut
    )
    if (is.null(ego) || nrow(ego) == 0) return(NULL)
    ego
}

enrich_all <- function(gene_list, direction) {
    res <- list()
    for (ont in ont_list) {
        ego <- run_enrich_go(gene_list, ont)
        if (!is.null(ego)) res[[ont]] <- ego
    }
    res
}

combine_ego_to_df <- function(ego_list, direction) {
    if (length(ego_list) == 0) return(NULL)
    df_list <- list()
    for (ont in names(ego_list)) {
        df <- as.data.frame(ego_list[[ont]])
        df$Direction <- direction
        df$Ontology <- ont
        df_list[[ont]] <- df
    }
    do.call(rbind, df_list)
}

deg <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)
required_cols <- c("GeneID", "log2FoldChange", "pvalue", "padj")
missing_cols <- setdiff(required_cols, colnames(deg))
if (length(missing_cols) > 0) {
    stop("DEG result is missing columns: ", paste(missing_cols, collapse = ", "))
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

for (sig_col in sig_rounds) {
    message("filter: ", sig_col, " < ", sig_threshold)
    sig_mask <- deg[, sig_col] < sig_threshold & !is.na(deg[, sig_col])
    fc_up <- deg$GeneID[sig_mask & deg$log2FoldChange > 0]
    fc_down <- deg$GeneID[sig_mask & deg$log2FoldChange < 0]

    up_df <- combine_ego_to_df(enrich_all(fc_up, "Up"), "Up")
    down_df <- combine_ego_to_df(enrich_all(fc_down, "Down"), "Down")

    suffix <- ifelse(sig_col == "padj", "_padj", "_pvalue")
    csv_file <- file.path(output_dir,
                          paste0(comp_name, "_GO_enrichment", suffix, ".csv"))
    xlsx_file <- file.path(output_dir,
                           paste0(comp_name, "_GO_enrichment", suffix, ".xlsx"))

    all_df <- rbind(up_df, down_df)
    if (!is.null(all_df) && nrow(all_df) > 0) {
        all_df <- all_df[order(all_df$Direction, all_df$Ontology,
                               all_df$p.adjust), ]
        all_df <- all_df[, c("Direction", "Ontology",
                             setdiff(colnames(all_df),
                                     c("Direction", "Ontology")))]
        write.csv(all_df, file = csv_file, row.names = FALSE)
        message("  wrote ", csv_file)
    } else {
        all_df <- NULL
    }

    wb <- createWorkbook()
    for (sheet in c("Up_regulated", "Down_regulated")) {
        addWorksheet(wb, sheet)
        df <- if (sheet == "Up_regulated") up_df else down_df
        if (!is.null(df) && nrow(df) > 0) {
            writeData(wb, sheet, df)
        } else {
            writeData(wb, sheet, "no enriched term")
        }
    }
    saveWorkbook(wb, xlsx_file, overwrite = TRUE)
    message("  wrote ", xlsx_file)
}
