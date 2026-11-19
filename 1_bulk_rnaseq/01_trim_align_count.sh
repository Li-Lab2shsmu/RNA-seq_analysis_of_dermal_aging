#!/bin/bash
#
# Purpose: adapter trimming, alignment and gene-level counting for the bulk
#          RNA-seq samples.
# Input:   data/fastq/<sample>.R1.fq.gz and <sample>.R2.fq.gz
# Output:  data/trimmed/<sample>.R1_val_1.fq.gz, <sample>.R2_val_2.fq.gz
#          data/mapping/<sample>.sorted.bam, <sample>.unmapped.bam,
#                       <sample>.hisat2_summary.txt
#          data/counts/all_samples.featureCounts.txt and .summary.txt
# Usage:   bash 01_trim_align_count.sh [trim|align|count|all]   (default: all)
#
# Each stage only needs the previous stage's output on disk, so a run that
# stopped halfway can be resumed from where it did. trim_galore, hisat2,
# samtools and featureCounts must be in PATH. Per-sample failures are reported
# at the end and make the run exit non-zero.

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# --- settings ---------------------------------------------------------------

# All sequenced samples, including D48. D48 is only removed downstream after QC.
SAMPLE_LIST=(D5 D8 D10 D21 D26 D30 D43 D48 D50 D62 D70 D74 D79)

DIR_FASTQ="data/fastq"
DIR_TRIMMED="data/trimmed"
DIR_MAPPING="data/mapping"
DIR_COUNTS="data/counts"

RAW_R1_SUFFIX=".R1.fq.gz"
RAW_R2_SUFFIX=".R2.fq.gz"
TRIMMED_R1_SUFFIX=".R1_val_1.fq.gz"
TRIMMED_R2_SUFFIX=".R2_val_2.fq.gz"

TRIM_MIN_LENGTH=20
TRIM_CLIP_5=10
TRIM_CLIP_3=10
THREADS=10

# Reference files differ between machines. Set them in the environment or edit
# the defaults before running.
HISAT2_INDEX="${HISAT2_INDEX:-/path/to/hg19/genome}"
GTF="${GTF:-/path/to/genome.gtf}"

LOG_FILE=""

log() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOG_FILE"
}

die() {
    printf '[%s] ERROR: %s\n' "$(date '+%F %T')" "$*" >&2
    exit 1
}

require_tools() {
    local tool
    for tool in "$@"; do
        command -v "$tool" >/dev/null 2>&1 || die "required tool not in PATH: $tool"
    done
}

require_tools trim_galore hisat2 samtools featureCounts

for ref in "$HISAT2_INDEX" "$GTF"; do
    [[ "$ref" == /path/to/* ]] && die "set HISAT2_INDEX and GTF before running"
done

failures=0


# --- trimming ---------------------------------------------------------------

run_trim() {
    LOG_FILE="${DIR_TRIMMED}/trim_galore.log"
    mkdir -p "$DIR_TRIMMED"
    log "trim_galore: ${#SAMPLE_LIST[@]} samples -> $DIR_TRIMMED" \
        "(min length $TRIM_MIN_LENGTH, clip ${TRIM_CLIP_5}/${TRIM_CLIP_3} bp)"

    local sample r1 r2 out_r1 out_r2
    for sample in "${SAMPLE_LIST[@]}"; do
        r1="${DIR_FASTQ}/${sample}${RAW_R1_SUFFIX}"
        r2="${DIR_FASTQ}/${sample}${RAW_R2_SUFFIX}"
        out_r1="${DIR_TRIMMED}/${sample}${TRIMMED_R1_SUFFIX}"
        out_r2="${DIR_TRIMMED}/${sample}${TRIMMED_R2_SUFFIX}"

        if [[ ! -f "$r1" || ! -f "$r2" ]]; then
            log "$sample: input FASTQ missing, skipped"
            failures=$((failures + 1))
            continue
        fi

        log "trimming $sample"
        if ! trim_galore \
                --paired "$r1" "$r2" \
                --length "$TRIM_MIN_LENGTH" \
                --clip_R1 "$TRIM_CLIP_5" \
                --clip_R2 "$TRIM_CLIP_5" \
                --three_prime_clip_R1 "$TRIM_CLIP_3" \
                --three_prime_clip_R2 "$TRIM_CLIP_3" \
                --output_dir "$DIR_TRIMMED" \
                2>>"$LOG_FILE"; then
            log "$sample: trim_galore failed"
            failures=$((failures + 1))
            continue
        fi

        if [[ ! -s "$out_r1" || ! -s "$out_r2" ]]; then
            log "$sample: trim_galore produced no output"
            failures=$((failures + 1))
        fi
    done
}


# --- alignment --------------------------------------------------------------

run_align() {
    LOG_FILE="${DIR_MAPPING}/hisat2_mapping.log"
    mkdir -p "$DIR_MAPPING"
    log "hisat2: ${#SAMPLE_LIST[@]} samples -> $DIR_MAPPING" \
        "(index $HISAT2_INDEX, $THREADS threads)"

    local sample r1 r2 prefix
    for sample in "${SAMPLE_LIST[@]}"; do
        r1="${DIR_TRIMMED}/${sample}${TRIMMED_R1_SUFFIX}"
        r2="${DIR_TRIMMED}/${sample}${TRIMMED_R2_SUFFIX}"
        prefix="${DIR_MAPPING}/${sample}"

        if [[ ! -f "$r1" || ! -f "$r2" ]]; then
            log "$sample: trimmed FASTQ missing, skipped"
            failures=$((failures + 1))
            continue
        fi

        log "aligning $sample"
        if ! hisat2 \
                -x "$HISAT2_INDEX" \
                -p "$THREADS" \
                -1 "$r1" \
                -2 "$r2" \
                -S "${prefix}.sam" \
                --summary-file "${prefix}.hisat2_summary.txt" \
                --very-sensitive \
                2>>"$LOG_FILE" || [[ ! -s "${prefix}.sam" ]]; then
            log "$sample: hisat2 failed"
            failures=$((failures + 1))
            continue
        fi

        if ! samtools view -b "${prefix}.sam" > "${prefix}.bam" 2>>"$LOG_FILE" \
                || [[ ! -s "${prefix}.bam" ]]; then
            log "$sample: samtools view failed"
            failures=$((failures + 1))
            continue
        fi
        rm -f "${prefix}.sam"

        if ! samtools sort "${prefix}.bam" -o "${prefix}.sorted.bam" 2>>"$LOG_FILE" \
                || [[ ! -s "${prefix}.sorted.bam" ]]; then
            log "$sample: samtools sort failed"
            failures=$((failures + 1))
            continue
        fi
        rm -f "${prefix}.bam"

        # Flag 12 = read unmapped (4) + mate unmapped (8). An empty result
        # means everything mapped, which is not an error.
        if ! samtools view -b -f 12 "${prefix}.sorted.bam" > "${prefix}.unmapped.bam" \
                2>>"$LOG_FILE"; then
            log "$sample: extracting unmapped reads failed"
        fi
        log "$sample: done"
    done
}


# --- counting ---------------------------------------------------------------

run_count() {
    LOG_FILE="${DIR_COUNTS}/featureCounts_merged.log"
    mkdir -p "$DIR_COUNTS"

    local merged="${DIR_COUNTS}/all_samples.featureCounts.txt"
    local summary="${DIR_COUNTS}/all_samples.featureCounts.summary.txt"

    # One featureCounts call covers every sample, so a missing BAM would quietly
    # produce a matrix with fewer columns. Wait until all samples are on disk.
    local bam_list=() missing=() sample
    for sample in "${SAMPLE_LIST[@]}"; do
        if [[ -s "${DIR_MAPPING}/${sample}.sorted.bam" ]]; then
            bam_list+=("${DIR_MAPPING}/${sample}.sorted.bam")
        else
            missing+=("$sample")
        fi
    done
    if (( ${#missing[@]} > 0 )); then
        log "missing sorted BAM: ${missing[*]}"
        die "run the align stage first (BAMs are expected in $DIR_MAPPING)"
    fi

    log "featureCounts: ${#bam_list[@]} BAMs -> $merged"
    # -p paired-end, -t exon -g gene_name to count by gene symbol.
    if ! featureCounts \
            -T "$THREADS" \
            -p \
            -t exon \
            -g gene_name \
            -a "$GTF" \
            -o "$merged" \
            "${bam_list[@]}" \
            2>&1 | tee -a "$LOG_FILE"; then
        die "featureCounts failed, see $LOG_FILE"
    fi
    [[ -s "$merged" ]] || die "featureCounts wrote no counts to $merged"

    if [[ -f "${merged}.summary" ]]; then
        mv -f "${merged}.summary" "$summary"
        log "per-sample assignment summary: $summary"
    else
        log "no .summary file next to the counts table"
    fi

    log "counted $(($(wc -l < "$merged") - 1)) genes in ${#bam_list[@]} samples"
}


# --- entry point ------------------------------------------------------------

stage="${1:-all}"
case "$stage" in
    trim)  run_trim ;;
    align) run_align ;;
    count) run_count ;;
    all)   run_trim
           run_align
           run_count ;;
    *)     die "unknown stage: $stage (use trim, align, count or all)" ;;
esac

if (( failures > 0 )); then
    printf '[%s] %d sample(s) did not finish, check the stage logs\n' \
        "$(date '+%F %T')" "$failures" >&2
    exit 1
fi
printf '[%s] %s: finished\n' "$(date '+%F %T')" "$stage"
