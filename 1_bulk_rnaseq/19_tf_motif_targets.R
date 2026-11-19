#!/usr/bin/env Rscript

# Purpose: infer NR4A1, ETV4 and MYC target genes from RcisTarget motif
#          rankings.
# Input:   data/cistarget/hg19-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather
#          data/cistarget/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
#          results/MAT_RPKM_noD48.tsv
# Output:  results/motif_targets/<TF>_motif_targets.csv
#          results/motif_targets/<TF>_motif_summary.csv
#          results/motif_targets/all_motif_targets.csv
# Usage:   Rscript 19_tf_motif_targets.R
#
# These three TFs have no usable fibroblast ChIP-seq, so targets are called
# from sequence motifs instead. A gene is a candidate when its best rank for a
# directly annotated motif falls in the top rank_top_pct of the genome and at
# least min_motif_support motifs support it. Rank 0 means the database has no
# value for that gene-motif pair and is excluded, not treated as rank 0.

suppressPackageStartupMessages(library(arrow))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
    stop("cannot locate this script, run it as: Rscript 19_tf_motif_targets.R")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
setwd(script_dir)
source("config.R")

cistarget_dir <- file.path("data", "cistarget")
feather_file <- file.path(cistarget_dir,
    "hg19-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather")
motif_annot_file <- file.path(cistarget_dir,
    "motifs-v9-nr.hgnc-m0.001-o0.0.tbl")
rpkm_file <- MAT_RPKM
output_dir <- file.path(DIR_RESULTS, "motif_targets")

tfs <- c("NR4A1", "ETV4", "MYC")
rank_top_pct <- 0.03
direct_annot_only <- TRUE
exclude_sources <- c("tfdimers", "taipale_tf_pairs")
min_expr <- 1
min_motif_support <- 2
overlap_with_coexpr <- TRUE
coexpr_rho <- 0.6
chunk_size <- 2000

if (!file.exists(feather_file)) {
    stop("motif ranking database not found: ", feather_file,
         "\n  download from https://resources.aertslab.org/cistarget/")
}
if (!file.exists(motif_annot_file)) {
    stop("motif annotation table not found: ", motif_annot_file,
         "\n  download from https://resources.aertslab.org/cistarget/")
}

all_cols <- names(arrow::open_dataset(feather_file, format = "feather")$schema)
if (!"motifs" %in% all_cols) {
    stop("the feather file has no 'motifs' column")
}

motif_col <- arrow::read_feather(feather_file, col_select = "motifs")
db_motifs <- as.character(motif_col$motifs)
rm(motif_col)
gc()

gene_cols <- setdiff(all_cols, "motifs")
n_gene <- length(gene_cols)
rank_cut <- max(1L, round(rank_top_pct * n_gene))
message(length(db_motifs), " motifs x ", n_gene, " genes, rank cut ", rank_cut)

annot <- read.delim(motif_annot_file, header = TRUE, sep = "\t",
                    comment.char = "", stringsAsFactors = FALSE, quote = "")
colnames(annot)[1] <- "motif"

need_col <- c("motif", "gene_name", "description")
miss_col <- setdiff(need_col, colnames(annot))
if (length(miss_col) > 0) {
    stop("annotation table is missing columns: ",
         paste(miss_col, collapse = ", "))
}

annot_tf <- annot[annot$gene_name %in% tfs, , drop = FALSE]
if (direct_annot_only) {
    annot_tf <- annot_tf[grepl("directly annotated", annot_tf$description), ,
                         drop = FALSE]
}
annot_tf$source <- sub("__.*$", "", annot_tf$motif)
if (length(exclude_sources) > 0) {
    annot_tf <- annot_tf[!annot_tf$source %in% exclude_sources, , drop = FALSE]
}
annot_tf <- annot_tf[annot_tf$motif %in% db_motifs, , drop = FALSE]
if (nrow(annot_tf) == 0) {
    stop("no usable motif for the requested TFs")
}

idx_motif <- match(annot_tf$motif, db_motifs)
annot_tf <- annot_tf[!is.na(idx_motif), , drop = FALSE]
idx_motif <- idx_motif[!is.na(idx_motif)]

hit_list <- list()
starts <- seq(1, n_gene, by = chunk_size)

for (st in starts) {
    cols <- gene_cols[st:min(st + chunk_size - 1, n_gene)]
    ch <- arrow::read_feather(feather_file, col_select = cols)
    cm <- as.matrix(ch)
    rm(ch)

    sub <- cm[idx_motif, , drop = FALSE]
    rm(cm)

    hit_chunk <- do.call(rbind, lapply(seq_along(idx_motif), function(i) {
        v <- sub[i, ]
        # Rank 0 is a missing value in the database, not the top rank.
        k <- which(!is.na(v) & v >= 1L & v <= rank_cut)
        if (length(k) == 0) return(NULL)
        data.frame(TF = annot_tf$gene_name[i],
                   motif = annot_tf$motif[i],
                   gene = cols[k],
                   rank = as.integer(v[k]),
                   stringsAsFactors = FALSE)
    }))
    rm(sub)

    if (!is.null(hit_chunk)) hit_list[[length(hit_list) + 1]] <- hit_chunk
    rm(hit_chunk)
    gc()
    message("columns ", st, "-", min(st + chunk_size - 1, n_gene), " done")
}

if (length(hit_list) == 0) {
    stop("no gene passed the rank cutoff, lower rank_top_pct")
}

hits <- do.call(rbind, hit_list)
rm(hit_list)
gc()

best_rank <- aggregate(rank ~ TF + gene, data = hits, FUN = min)
names(best_rank)[3] <- "best_rank"

mot_str <- aggregate(motif ~ TF + gene, data = hits,
                     FUN = function(x) paste(sort(unique(x)), collapse = ";"))
n_mot <- aggregate(motif ~ TF + gene, data = hits,
                   FUN = function(x) length(unique(x)))
names(n_mot)[3] <- "n_motifs"

res_all <- merge(best_rank, mot_str, by = c("TF", "gene"))
res_all <- merge(res_all, n_mot, by = c("TF", "gene"))
names(res_all)[names(res_all) == "motif"] <- "motifs"
res_all$rank_pct <- round(res_all$best_rank / n_gene, 5)

n_motif_of <- table(annot_tf$gene_name)

if (!file.exists(rpkm_file)) {
    stop("RPKM matrix not found: ", rpkm_file)
}
rpkm <- read.table(rpkm_file, header = TRUE, row.names = 1, sep = "\t",
                   check.names = FALSE, stringsAsFactors = FALSE)
mean_expr <- rowMeans(as.matrix(rpkm), na.rm = TRUE)

res_all$mean_expr <- round(mean_expr[res_all$gene], 3)
res_all$expressed <- res_all$gene %in% names(mean_expr)[mean_expr > min_expr]
res_all <- res_all[res_all$expressed, , drop = FALSE]

if (overlap_with_coexpr) {
    res_all$rho_vs_TF <- NA_real_
    for (tf in unique(res_all$TF)) {
        if (!tf %in% rownames(rpkm)) next
        tf_vec <- as.numeric(rpkm[tf, ])
        idx <- which(res_all$TF == tf)
        res_all$rho_vs_TF[idx] <- apply(
            rpkm[res_all$gene[idx], , drop = FALSE], 1,
            function(x) suppressWarnings(
                cor(as.numeric(x), tf_vec, method = "spearman")))
    }
    res_all$rho_vs_TF <- round(res_all$rho_vs_TF, 3)
    res_all$coexpressed <- !is.na(res_all$rho_vs_TF) &
        res_all$rho_vs_TF > coexpr_rho
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
all_out <- list()

for (tf in tfs) {
    res <- res_all[res_all$TF == tf, , drop = FALSE]
    if (nrow(res) == 0) next

    need_n <- min(min_motif_support, as.integer(n_motif_of[tf]))
    res <- res[res$n_motifs >= need_n, , drop = FALSE]
    res <- res[order(res$best_rank), , drop = FALSE]
    if (nrow(res) == 0) next
    res <- res[, c("TF", "gene", "best_rank", "rank_pct", "n_motifs", "motifs",
                   "mean_expr", "expressed",
                   intersect(c("rho_vs_TF", "coexpressed"), colnames(res))),
               drop = FALSE]

    write.csv(res,
              file = file.path(output_dir, paste0(tf, "_motif_targets.csv")),
              row.names = FALSE)

    all_out[[tf]] <- res

    mot_df <- annot_tf[annot_tf$gene_name == tf, c("motif", "source")]
    mot_df$n_targets <- sapply(mot_df$motif,
        function(m) sum(grepl(m, res$motifs, fixed = TRUE)))
    write.csv(mot_df,
              file = file.path(output_dir, paste0(tf, "_motif_summary.csv")),
              row.names = FALSE)
    message(tf, ": ", nrow(res), " candidate targets")
}

if (length(all_out) == 0) {
    stop("no TF passed the filters, lower rank_top_pct or min_motif_support")
}

combined <- do.call(rbind, all_out)
common_col <- Reduce(intersect, lapply(all_out, colnames))
write.csv(combined[, common_col, drop = FALSE],
          file = file.path(output_dir, "all_motif_targets.csv"),
          row.names = FALSE)
message("wrote ", file.path(output_dir, "all_motif_targets.csv"))
