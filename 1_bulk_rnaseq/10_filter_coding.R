#!/usr/bin/env Rscript

# Purpose: keep only protein-coding genes in the count matrix.
# Input:   results/MAT_COUNT_noD48.tsv, a GTF annotation
# Output:  results/coding_filter/protein_coding_genes.txt
#          results/coding_filter/MAT_COUNT_coding_noD48.tsv
# Usage:   GTF=/path/to/genome.gtf Rscript 10_filter_coding.R
#
# A gene name is kept when any GTF record with that name is protein_coding,
# so genes with incomplete transcript annotation are not lost.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
    stop("cannot locate this script, run it as: Rscript 10_filter_coding.R")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
setwd(script_dir)
source("config.R")

gtf_file <- Sys.getenv("GTF", unset = "/path/to/genome.gtf")
counts_file <- MAT_COUNT
output_dir <- file.path(DIR_RESULTS, "coding_filter")

if (gtf_file == "/path/to/genome.gtf") {
    stop("set GTF to the annotation file before running")
}
if (!file.exists(gtf_file)) {
    stop("GTF not found: ", gtf_file)
}
if (!file.exists(counts_file)) {
    stop("count matrix not found: ", counts_file,
         "\n  run 03_expression_matrix.R first")
}

message("reading GTF attribute column")
gtf_col9 <- read.table(gtf_file, header = FALSE, sep = "\t",
                       quote = "", comment.char = "#",
                       stringsAsFactors = FALSE)[, 9]
message(length(gtf_col9), " GTF records")

name_m <- regexpr('gene_name "([^"]+)"', gtf_col9, perl = TRUE)
gene_name <- ifelse(name_m > 0,
                    gsub('gene_name "|"$', "", regmatches(gtf_col9, name_m)),
                    NA_character_)

bio_m <- regexpr('gene_biotype "([^"]+)"', gtf_col9, perl = TRUE)
gene_biotype <- ifelse(bio_m > 0,
                       gsub('gene_biotype "|"$', "",
                            regmatches(gtf_col9, bio_m)),
                       NA_character_)

coding_mask <- gene_biotype == "protein_coding" & !is.na(gene_name)
coding_genes <- unique(gene_name[coding_mask])
coding_genes <- coding_genes[!is.na(coding_genes)]
message(length(coding_genes), " protein-coding gene names")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.table(coding_genes,
            file = file.path(output_dir, "protein_coding_genes.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

counts <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t",
                     stringsAsFactors = FALSE, check.names = FALSE)
keep_genes <- intersect(rownames(counts), coding_genes)
counts_coding <- counts[keep_genes, , drop = FALSE]
message(nrow(counts), " genes -> ", nrow(counts_coding), " protein-coding genes")

write.table(data.frame(GeneID = rownames(counts_coding), counts_coding,
                       check.names = FALSE),
            file = file.path(output_dir, "MAT_COUNT_coding_noD48.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote ", file.path(output_dir, "MAT_COUNT_coding_noD48.tsv"))
