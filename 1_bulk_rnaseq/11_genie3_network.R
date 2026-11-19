#!/usr/bin/env Rscript

# Purpose: infer TF-target regulatory networks with GENIE3 for the Middle and
#          Aged groups.
# Input:   results/MAT_RPKM_noD48.tsv
#          a protein-coding GTF set through the PC_GTF environment variable
# Output:  results/genie3/GENIE3_linklist_{Middle,Aged}.tsv
# Usage:   PC_GTF=/path/to/protein_coding.gtf Rscript 11_genie3_network.R
#
# Genes are log2(RPKM + 1) transformed. Targets are the most variable genes
# within each group, regulators are annotated transcription factors. With three
# samples per group the networks are exploratory.

suppressPackageStartupMessages({
    library(GENIE3)
    library(org.Hs.eg.db)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
    stop("cannot locate this script, run it as: Rscript 11_genie3_network.R")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
setwd(script_dir)
source("config.R")

rpkm_file <- MAT_RPKM
output_dir <- file.path(DIR_RESULTS, "genie3")
pc_gtf <- Sys.getenv("PC_GTF", unset = "")

n_top_var <- 3000
n_trees <- 1000
n_cores <- as.integer(Sys.getenv("GENIE3_CORES", unset = "8"))

if (!file.exists(rpkm_file)) {
    stop("RPKM matrix not found: ", rpkm_file,
         "\n  run 03_expression_matrix.R first")
}
if (pc_gtf == "" || !file.exists(pc_gtf)) {
    stop("set PC_GTF to a protein-coding GTF file before running")
}

rpkm <- read.table(rpkm_file, header = TRUE, row.names = 1, sep = "\t",
                   check.names = FALSE)
expr_all <- log2(rpkm + 1)
message(nrow(expr_all), " genes x ", ncol(expr_all), " samples")

pc_names <- system(
    sprintf("grep -oP 'gene_name \"[^\"]+\"' '%s'", pc_gtf),
    intern = TRUE
)
pc_names <- unique(sub('gene_name "', "", pc_names, fixed = TRUE))
pc_names <- sub('"$', "", pc_names)
expr_all <- expr_all[rownames(expr_all) %in% pc_names, , drop = FALSE]
message(nrow(expr_all), " protein-coding genes")

tf_ids <- unique(c(
    AnnotationDbi::select(org.Hs.eg.db, keys = "GO:0003700",
                          columns = "SYMBOL", keytype = "GOALL")$SYMBOL,
    AnnotationDbi::select(org.Hs.eg.db, keys = "GO:0000981",
                          columns = "SYMBOL", keytype = "GOALL")$SYMBOL
))
tf_set <- unique(tf_ids[!is.na(tf_ids)])
message(length(tf_set), " annotated transcription factors")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

for (net in GENIE3_NETWORKS) {
    samples <- GROUPS[[net]]
    missing <- setdiff(samples, colnames(expr_all))
    if (length(missing) > 0) {
        stop("samples missing from the matrix: ",
             paste(missing, collapse = ", "))
    }

    expr <- as.matrix(expr_all[, samples, drop = FALSE])

    # Drop genes that stay near the detection floor in every sample.
    keep <- rowSums(expr > 1) >= 1
    expr <- expr[keep, , drop = FALSE]

    gene_var <- apply(expr, 1, var)
    targets <- names(sort(gene_var, decreasing = TRUE))
    targets <- targets[seq_len(min(n_top_var, length(targets)))]
    regulators <- intersect(tf_set, rownames(expr))

    message(net, ": ", length(targets), " targets, ",
            length(regulators), " regulators")

    set.seed(42)
    weight_mat <- GENIE3(
        expr,
        regulators = regulators,
        targets    = targets,
        nTrees     = n_trees,
        nCores     = n_cores,
        verbose    = TRUE
    )

    link_list <- getLinkList(weight_mat)
    out_file <- file.path(output_dir, paste0("GENIE3_linklist_", net, ".tsv"))
    write.table(link_list, file = out_file,
                sep = "\t", quote = FALSE, row.names = FALSE)
    message("  wrote ", out_file, " (", nrow(link_list), " links)")
}
