# Shared settings for the bulk RNA-seq pipeline.
# Source from each script: source(file.path(script_dir, "config.R"))

# Sequencing samples, including D48; D48 is dropped downstream after QC.
ALL_SAMPLES <- c("D5", "D8", "D10", "D21", "D26", "D30", "D43",
                 "D48", "D50", "D62", "D70", "D74", "D79")
DROP_SAMPLES <- c("D48")

# Age groups used by Mfuzz and DESeq2 (D48 already removed).
GROUPS <- list(
    Young  = c("D5", "D8", "D10"),
    Adult  = c("D21", "D26", "D30"),
    Middle = c("D43", "D50", "D62"),
    Aged   = c("D70", "D74", "D79")
)

# Differential expression: Aged is the test group, Middle the control.
DEG_TEST <- "Aged"
DEG_CONTROL <- "Middle"

# Networks inferred by GENIE3.
GENIE3_NETWORKS <- c("Middle", "Aged")

# Directories, relative to 1_bulk_rnaseq/.
DIR_FASTQ <- "data/fastq"
DIR_TRIMMED <- "data/trimmed"
DIR_MAPPING <- "data/mapping"
DIR_COUNTS <- "data/counts"
DIR_QC <- "results/qc"
DIR_RESULTS <- "results"

# Count and RPKM matrices read by the downstream scripts.
MAT_COUNT <- file.path(DIR_RESULTS, "MAT_COUNT_noD48.tsv")
MAT_RPKM <- file.path(DIR_RESULTS, "MAT_RPKM_noD48.tsv")

# Significance cutoff used for DEG tables and downstream screens.
SIG_CRITERIA <- "pvalue"
SIG_THRESHOLD <- 0.05
MIN_COUNT <- 10
