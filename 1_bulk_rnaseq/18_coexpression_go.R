#!/usr/bin/env Rscript

# Purpose: GO enrichment of genes co-expressed with NR4A1, MYC and ETV4.
# Input:   results/MAT_RPKM_noD48.tsv
#          a protein-coding GTF set through the PC_GTF environment variable
# Output:  results/genie3_go/<TF>_coexpression_GO.csv
#          results/genie3_go/all_coexpression_GO.csv
# Usage:   PC_GTF=/path/to/protein_coding.gtf Rscript 18_coexpression_go.R
#
# Co-expressed genes are protein-coding genes with Spearman rho > 0.6 against
# the TF across the 12 samples, which captures programmes that move with it.

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(dplyr)
})

options(warn = -1)
on.exit(options(warn = 0))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
    stop("cannot locate this script, run it as: Rscript 18_coexpression_go.R")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
setwd(script_dir)
source("config.R")

rpkm_file <- MAT_RPKM
pc_gtf <- Sys.getenv("PC_GTF", unset = "")
output_dir <- file.path(DIR_RESULTS, "genie3_go")

tfs <- c("NR4A1", "MYC", "ETV4")
rho_threshold <- 0.6
min_expr <- 1

species_db <- "org.Hs.eg.db"
key_type <- "SYMBOL"
ont <- "BP"
pvalue_cut <- 0.05
qvalue_cut <- 0.2

if (!file.exists(rpkm_file)) {
    stop("RPKM matrix not found: ", rpkm_file,
         "\n  run 03_expression_matrix.R first")
}
if (pc_gtf == "" || !file.exists(pc_gtf)) {
    stop("set PC_GTF to a protein-coding GTF file before running")
}

rpkm <- read.table(rpkm_file, header = TRUE, row.names = 1, sep = "\t",
                   check.names = FALSE)
rpkm <- rpkm[rowMeans(rpkm) > min_expr, , drop = FALSE]

pc_names <- system(sprintf("grep -oP 'gene_name \"[^\"]+\"' '%s'", pc_gtf),
                   intern = TRUE)
pc_names <- unique(sub('gene_name "', "", pc_names, fixed = TRUE))
pc_names <- sub('"$', "", pc_names)
rpkm <- rpkm[rownames(rpkm) %in% pc_names, , drop = FALSE]
message(nrow(rpkm), " protein-coding genes above the expression filter")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
all_results <- list()

for (tf in tfs) {
    if (!tf %in% rownames(rpkm)) {
        message(tf, ": not in the expression matrix, skipped")
        next
    }

    tf_expr <- as.numeric(rpkm[tf, ])
    rho_vec <- apply(rpkm, 1, function(x) {
        suppressWarnings(cor(tf_expr, as.numeric(x), method = "spearman"))
    })

    genes <- names(rho_vec)[rho_vec > rho_threshold & !is.na(rho_vec)]
    genes <- setdiff(genes, tf)
    message(tf, ": ", length(genes), " co-expressed genes")

    if (length(genes) < 5) {
        message(tf, ": fewer than 5 genes, skipped")
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
              file = file.path(output_dir, paste0(tf, "_coexpression_GO.csv")),
              row.names = FALSE)
}

if (length(all_results) > 0) {
    combined <- do.call(rbind, all_results)
    write.csv(combined, file = file.path(output_dir, "all_coexpression_GO.csv"),
              row.names = FALSE)
    message("wrote ", file.path(output_dir, "all_coexpression_GO.csv"))
}
