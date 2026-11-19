#!/bin/bash
#
# Purpose: quality control of the trimmed reads and of the alignments,
#          collected into one MultiQC report.
# Input:   data/trimmed/<sample>.R1_val_1.fq.gz, <sample>.R2_val_2.fq.gz
#          data/mapping/<sample>.sorted.bam
# Output:  results/qc/fastqc/*_fastqc.html
#          results/qc/bam_stats/<sample>.stats.txt
#          results/qc/multiqc/combined_qc_report.html
# Usage:   bash 02_qc.sh
#
# fastqc, samtools and multiqc must be in PATH. QC is diagnostic, so a failed
# sample is reported and skipped; only a complete failure of both stages stops
# the run.

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# --- settings ---------------------------------------------------------------

SAMPLE_LIST=(D5 D8 D10 D21 D26 D30 D43 D48 D50 D62 D70 D74 D79)

DIR_TRIMMED="data/trimmed"
DIR_MAPPING="data/mapping"
DIR_QC="results/qc"

TRIMMED_R1_SUFFIX=".R1_val_1.fq.gz"
TRIMMED_R2_SUFFIX=".R2_val_2.fq.gz"
QC_THREADS=4

LOG_FILE="${DIR_QC}/qc_pipeline.log"

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

require_tools fastqc samtools multiqc

FASTQC_DIR="${DIR_QC}/fastqc"
BAMSTATS_DIR="${DIR_QC}/bam_stats"
MULTIQC_DIR="${DIR_QC}/multiqc"
mkdir -p "$FASTQC_DIR" "$BAMSTATS_DIR" "$MULTIQC_DIR"

fastqc_ok=0
stats_ok=0


# --- FastQC on the trimmed reads --------------------------------------------

run_fastqc() {
    log "FastQC: ${#SAMPLE_LIST[@]} samples -> $FASTQC_DIR"

    local sample r1 r2 html1 html2
    for sample in "${SAMPLE_LIST[@]}"; do
        r1="${DIR_TRIMMED}/${sample}${TRIMMED_R1_SUFFIX}"
        r2="${DIR_TRIMMED}/${sample}${TRIMMED_R2_SUFFIX}"

        if [[ ! -f "$r1" || ! -f "$r2" ]]; then
            log "$sample: trimmed FASTQ missing, skipped"
            continue
        fi

        # FastQC names its report after the input file without the .fq.gz part.
        html1="${FASTQC_DIR}/${sample}${TRIMMED_R1_SUFFIX%.fq.gz}_fastqc.html"
        html2="${FASTQC_DIR}/${sample}${TRIMMED_R2_SUFFIX%.fq.gz}_fastqc.html"

        if ! fastqc \
                --threads "$QC_THREADS" \
                --outdir "$FASTQC_DIR" \
                --quiet \
                "$r1" "$r2" \
                2>>"$LOG_FILE"; then
            log "$sample: fastqc failed"
            continue
        fi

        if [[ -f "$html1" && -f "$html2" ]]; then
            fastqc_ok=$((fastqc_ok + 1))
        else
            log "$sample: fastqc produced no report"
        fi
    done

    log "FastQC finished for $fastqc_ok/${#SAMPLE_LIST[@]} samples"
}


# --- samtools stats on the sorted BAMs --------------------------------------

run_bam_stats() {
    log "samtools stats: ${#SAMPLE_LIST[@]} samples -> $BAMSTATS_DIR"

    local sample bam out
    for sample in "${SAMPLE_LIST[@]}"; do
        bam="${DIR_MAPPING}/${sample}.sorted.bam"
        out="${BAMSTATS_DIR}/${sample}.stats.txt"

        if [[ ! -f "$bam" ]]; then
            log "$sample: sorted BAM missing, skipped"
            continue
        fi

        if ! samtools stats -@ "$QC_THREADS" "$bam" > "$out" 2>>"$LOG_FILE" \
                || [[ ! -s "$out" ]]; then
            log "$sample: samtools stats failed"
            continue
        fi
        stats_ok=$((stats_ok + 1))
    done

    log "samtools stats finished for $stats_ok/${#SAMPLE_LIST[@]} samples"
}


# --- MultiQC over both result directories -----------------------------------

run_multiqc() {
    local report="${MULTIQC_DIR}/combined_qc_report.html"
    log "MultiQC: scanning $FASTQC_DIR and $BAMSTATS_DIR"

    if ! multiqc \
            "$FASTQC_DIR" \
            "$BAMSTATS_DIR" \
            --outdir "$MULTIQC_DIR" \
            --filename "combined_qc_report" \
            --title "RNA-seq Quality Report (FASTQ + BAM)" \
            2>&1 | tee -a "$LOG_FILE"; then
        die "MultiQC failed, see $LOG_FILE"
    fi
    [[ -f "$report" ]] || die "MultiQC report not found: $report"

    log "report: $report"
}


# --- entry point ------------------------------------------------------------

run_fastqc
run_bam_stats

if (( fastqc_ok == 0 && stats_ok == 0 )); then
    die "both FastQC and samtools stats failed, nothing to summarise"
fi

run_multiqc
printf '[%s] qc: finished\n' "$(date '+%F %T')"
