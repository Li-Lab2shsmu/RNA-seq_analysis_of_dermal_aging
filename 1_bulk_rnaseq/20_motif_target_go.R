#!/usr/bin/env Rscript

# Purpose: GO enrichment of the motif-inferred targets of NR4A1, ETV4 and MYC.
# Input:   results/motif_targets/all_motif_targets.csv
#          a protein-coding GTF set through the PC_GTF environment variable
# Output:  results/motif_targets/target_go/<TF>_motif_GO.csv
#          results/motif_targets/target_go/<TF>_motif_coexpr_GO.csv
#          results/motif_targets/target_go/all_motif_GO.csv
# Usage:   PC_GTF=/path/to/protein_coding.gtf Rscript 20_motif_target_go.R
#
# Two gene sets are tested per TF: all motif candidates, and only candidates
# co-expressed with the TF. The second set combines sequence and expression
# evidence and is the stricter one.

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
})

options(warn = -1)
on.exit(options(warn = 0))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
    stop("cannot locate this script, run it as: Rscript 20_motif_target_go.R")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
setwd(script_dir)
source("config.R")

input_dir <- file.path(DIR_RESULTS, "motif_targets")
input_file <- file.path(input_dir, "all_motif_targets.csv")
pc_gtf <- Sys.getenv("PC_GTF", unset = "")
output_dir <- file.path(input_dir, "target_go")

tfs <- c("NR4A1", "ETV4", "MYC")
run_all_targets <- TRUE
run_coexpr_only <- TRUE
min_gene <- 5

species_db <- "org.Hs.eg.db"
key_type <- "SYMBOL"
ont <- "BP"
pvalue_cut <- 0.05
qvalue_cut <- 0.2

if (!file.exists(input_file)) {
    stop("motif target table not found: ", input_file,
         "\n  run 19_tf_motif_targets.R first")
}
if (pc_gtf == "" || !file.exists(pc_gtf)) {
    stop("set PC_GTF to a protein-coding GTF file before running")
}

tg <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)
if (!all(c("TF", "gene") %in% colnames(tg))) {
    stop("motif target table is missing TF or gene columns")
}
if (run_coexpr_only && !"coexpressed" %in% colnames(tg)) {
    stop("no coexpressed column; set run_coexpr_only to FALSE")
}

pc_names <- system(sprintf("grep -oP 'gene_name \"[^\"]+\"' '%s'", pc_gtf),
                   intern = TRUE)
pc_names <- unique(sub('gene_name "', "", pc_names, fixed = TRUE))
pc_names <- sub('"$', "", pc_names)
tg <- tg[tg$gene %in% pc_names, , drop = FALSE]

known_symbols <- AnnotationDbi::keys(org.Hs.eg.db, keytype = key_type)

run_go <- function(genes, tag) {
    genes <- unique(genes)
    genes <- genes[genes %in% known_symbols]
    if (length(genes) < min_gene) {
        message(tag, ": only ", length(genes), " usable genes, skipped")
        return(NULL)
    }
    ego <- enrichGO(gene = genes, OrgDb = species_db, keyType = key_type,
                    ont = ont, pAdjustMethod = "BH",
                    pvalueCutoff = pvalue_cut, qvalueCutoff = qvalue_cut)
    if (is.null(ego) || nrow(ego) == 0) {
        message(tag, ": no enriched term")
        return(NULL)
    }
    res <- as.data.frame(ego)
    res <- res[order(res$p.adjust), , drop = FALSE]
    res$TF <- NA_character_
    res$geneSet <- tag
    res$n_input <- length(genes)
    res
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
all_results <- list()

for (tf in tfs) {
    sub <- tg[tg$TF == tf, , drop = FALSE]
    if (nrow(sub) == 0) next

    if (run_all_targets) {
        r <- run_go(sub$gene, "motif")
        if (!is.null(r)) {
            r$TF <- tf
            all_results[[paste(tf, "motif", sep = "_")]] <- r
            write.csv(r, file = file.path(output_dir,
                                          paste0(tf, "_motif_GO.csv")),
                      row.names = FALSE)
        }
    }

    if (run_coexpr_only) {
        r2 <- run_go(sub$gene[sub$coexpressed], "motif_coexpr")
        if (!is.null(r2)) {
            r2$TF <- tf
            all_results[[paste(tf, "coexpr", sep = "_")]] <- r2
            write.csv(r2, file = file.path(output_dir,
                                           paste0(tf, "_motif_coexpr_GO.csv")),
                      row.names = FALSE)
        }
    }
}

if (length(all_results) > 0) {
    combined <- do.call(rbind, all_results)
    combined <- combined[, c("TF", "geneSet",
                             colnames(combined)[!colnames(combined) %in%
                                                c("TF", "geneSet")])]
    write.csv(combined, file = file.path(output_dir, "all_motif_GO.csv"),
              row.names = FALSE)
    message("wrote ", file.path(output_dir, "all_motif_GO.csv"))
} else {
    message("no gene set produced an enriched term")
}
