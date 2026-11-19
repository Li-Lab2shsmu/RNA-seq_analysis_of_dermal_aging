#!/usr/bin/env Rscript

# Purpose: turn the featureCounts table into the count and RPKM matrices read
#          by the downstream scripts.
# Input:   data/counts/all_samples.featureCounts.txt
# Output:  results/MAT_COUNT.tsv, results/MAT_RPKM.tsv
#          results/MAT_COUNT_noD48.tsv, results/MAT_RPKM_noD48.tsv
# Usage:   Rscript 03_expression_matrix.R
#
# D48 is kept in the full matrices and removed from the *_noD48 matrices. The
# noD48 names are fixed because every downstream script reads them.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
    stop("cannot locate this script, run it as: Rscript 03_expression_matrix.R")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
setwd(script_dir)
source("config.R")

counts_file <- file.path(DIR_COUNTS, "all_samples.featureCounts.txt")
if (!file.exists(counts_file)) {
    stop("featureCounts table not found: ", counts_file,
         "\n  run 01_trim_align_count.sh first")
}
dir.create(DIR_RESULTS, showWarnings = FALSE, recursive = TRUE)


# "D5.sorted.bam" -> "D5". basename() first so the rule also holds if the BAM
# paths were passed in as absolute paths.
clean_sample_names <- function(x) {
    sub("\\.sorted\\.bam$", "", basename(x))
}


# The first line is featureCounts' command banner, the second is the header.
# After Geneid the columns are Chr, Start, End, Strand, Length and then one
# column per sample; dropping the four annotation columns leaves Length in
# front of the counts. row.names = 1 errors out on duplicated gene symbols,
# which is intended: silently renaming them would break every downstream join
# by gene name.
read_featurecounts <- function(path) {
    raw <- read.table(path, header = TRUE, skip = 1, sep = "\t",
                      stringsAsFactors = FALSE, row.names = 1,
                      check.names = FALSE)

    annotation_cols <- c("Chr", "Start", "End", "Strand")
    missing <- setdiff(c(annotation_cols, "Length"), colnames(raw))
    if (length(missing) > 0) {
        stop("unexpected columns in ", basename(path), ": missing ",
             paste(missing, collapse = ", "))
    }
    raw <- raw[, setdiff(colnames(raw), annotation_cols), drop = FALSE]

    gene_length <- raw[[1]]
    names(gene_length) <- rownames(raw)
    counts <- as.matrix(raw[, -1, drop = FALSE])
    colnames(counts) <- clean_sample_names(colnames(counts))

    list(gene_length = gene_length, counts = counts)
}


# Genes with zero counts in every sample carry no information and would be
# all-zero rows in both matrices.
drop_unexpressed <- function(counts, gene_length) {
    keep <- rowSums(counts) > 0
    list(gene_length = gene_length[keep], counts = counts[keep, , drop = FALSE])
}


# RPKM = reads / gene length in kb / library size in millions.
to_rpkm <- function(counts, gene_length) {
    rpkm <- apply(counts, 2, function(x) {
        total <- sum(x)
        if (total == 0) return(x)
        x / gene_length * 1000 / total * 1e6
    })
    rownames(rpkm) <- rownames(counts)
    rpkm
}


remove_samples <- function(mat, samples) {
    missing <- setdiff(samples, colnames(mat))
    if (length(missing) > 0) {
        stop("cannot drop samples that are not in the matrix: ",
             paste(missing, collapse = ", "))
    }
    mat[, setdiff(colnames(mat), samples), drop = FALSE]
}


write_matrix <- function(mat, path) {
    out <- data.frame(GeneID = rownames(mat), mat, check.names = FALSE)
    write.table(out, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
    message("  ", path, "  (", nrow(mat), " genes x ", ncol(mat), " samples)")
}


main <- function() {
    message("reading ", counts_file)
    raw <- read_featurecounts(counts_file)
    n_raw <- nrow(raw$counts)
    message("  ", n_raw, " genes x ", ncol(raw$counts), " samples")
    message("  sample names: ", paste(colnames(raw$counts), collapse = ", "))

    kept <- drop_unexpressed(raw$counts, raw$gene_length)
    message("kept ", nrow(kept$counts), " genes (dropped ",
            n_raw - nrow(kept$counts), " with zero counts in every sample)")

    rpkm <- to_rpkm(kept$counts, kept$gene_length)

    message("writing matrices to ", DIR_RESULTS)
    write_matrix(kept$counts, file.path(DIR_RESULTS, "MAT_COUNT.tsv"))
    write_matrix(rpkm, file.path(DIR_RESULTS, "MAT_RPKM.tsv"))

    if (length(DROP_SAMPLES) > 0) {
        message("writing the matrices without ",
                paste(DROP_SAMPLES, collapse = ", "))
        write_matrix(remove_samples(kept$counts, DROP_SAMPLES),
                     file.path(DIR_RESULTS, "MAT_COUNT_noD48.tsv"))
        write_matrix(remove_samples(rpkm, DROP_SAMPLES),
                     file.path(DIR_RESULTS, "MAT_RPKM_noD48.tsv"))
    }
}

main()
