#!/bin/bash
# ==============================================
# SortMeRNA Batch Script - FAST VERSION
# ==============================================
# Using SortMeRNA v4.3 optimized database (99.9% accuracy)
# Estimated runtime: ~3-5 min per sample with 10 threads
#
# Usage:
#   conda activate pcc11901_rnaseq
#   bash logs/run_sortmerna_fast.sh
#
# To resume after interruption, just re-run - it auto-skips completed samples.
# ==============================================

set -e  # Exit on error

# Threading configuration (M2 Max: 12 cores, 64GB RAM)
# With 140MB database, memory is not a bottleneck - use 10 threads (leaves 2 for OS)
THREADS=10

# Base directory (handles spaces in path)
BASE_DIR="/Users/felix/Library/CloudStorage/OneDrive-SharedLibraries-MacquarieUniversity/Australian Genome Foundry - AWS cloud infrastructure/11_Esther_Cyano_transcriptomics"

# SortMeRNA v4.3 optimized database (recommended)
REF_DB="$BASE_DIR/DB/sortmerna_v4.3_db/smr_v4.3_default_db.fasta"

# Progress bar function (tqdm-style)
progress_bar() {
    local current=$1
    local total=$2
    local sample=$3
    local status=$4
    local eta=$5

    local width=30
    local percent=$((current * 100 / total))
    local filled=$((current * width / total))
    local empty=$((width - filled))

    # Build the bar
    local bar=""
    for ((i=0; i<filled; i++)); do bar+="█"; done
    for ((i=0; i<empty; i++)); do bar+="░"; done

    # Print progress line
    printf "\r%3d%%|%s| %d/%d [%s] %s %s" "$percent" "$bar" "$current" "$total" "$sample" "$status" "$eta"
}

# Verify database exists
if [ ! -f "$REF_DB" ]; then
    echo "ERROR: Database not found at $REF_DB"
    exit 1
fi

# Create output directories
mkdir -p "$BASE_DIR/04_rRNA_filtered/rRNA"
mkdir -p "$BASE_DIR/04_rRNA_filtered/non_rRNA"

echo "=============================================="
echo "SortMeRNA Batch Processing"
echo "=============================================="
echo "Database: SortMeRNA v4.3 default (99.9% accuracy)"
echo "Threads:  $THREADS"
echo "Started:  $(date '+%Y-%m-%d %H:%M:%S')"
echo "=============================================="
echo ""

# Sample list
SAMPLES=(U1 U2 U3 U4 U5 U6 U7 U8 U9 U10 U11 U12 U13 U14 U15 U16 U17 U18 U19 U20 U21 U22 U23 U24 U25 U26 U27 U28 U29 U30 U31 U32 U33 U34 U35 U36 U37 U38 U39 U40 U41 U42 U43 U44 U45 U46 U47 U48 U49 U50 U51 U52 U53 U54 U55 U56 U57 U58 U59 U60)

TOTAL=${#SAMPLES[@]}
COUNT=0
PROCESSED=0
SKIPPED=0
START_TIME=$(date +%s)
AVG_TIME=0

for SAMPLE in "${SAMPLES[@]}"; do
    COUNT=$((COUNT + 1))

    # Check if already processed (look for output files)
    if [ -f "$BASE_DIR/04_rRNA_filtered/non_rRNA/${SAMPLE}_fwd.fq.gz" ] || \
       [ -f "$BASE_DIR/04_rRNA_filtered/non_rRNA/${SAMPLE}_fwd.fq" ]; then
        SKIPPED=$((SKIPPED + 1))
        progress_bar $COUNT $TOTAL "$SAMPLE" "skipped" ""
        echo ""
        continue
    fi

    # Show progress bar with "running" status
    if [ $AVG_TIME -gt 0 ]; then
        REMAINING=$((TOTAL - COUNT))
        ETA_MINS=$(( (AVG_TIME * REMAINING) / 60 ))
        progress_bar $COUNT $TOTAL "$SAMPLE" "running..." "ETA: ${ETA_MINS}m"
    else
        progress_bar $COUNT $TOTAL "$SAMPLE" "running..." ""
    fi

    SAMPLE_START=$(date +%s)

    # Create workdir for this sample
    WORKDIR="$BASE_DIR/04_rRNA_filtered/workdir_${SAMPLE}"
    mkdir -p "$WORKDIR"

    # Run SortMeRNA (suppress output, log to file)
    sortmerna \
        --ref "$REF_DB" \
        --reads "$BASE_DIR/03_merged/${SAMPLE}_R1.fastq.gz" \
        --reads "$BASE_DIR/03_merged/${SAMPLE}_R2.fastq.gz" \
        --paired_in --out2 \
        --aligned "$BASE_DIR/04_rRNA_filtered/rRNA/${SAMPLE}" \
        --other "$BASE_DIR/04_rRNA_filtered/non_rRNA/${SAMPLE}" \
        --fastx \
        --threads $THREADS \
        --workdir "$WORKDIR" \
        > "$BASE_DIR/logs/${SAMPLE}_sortmerna.log" 2>&1

    # Clean up workdir (contains large index files)
    rm -rf "$WORKDIR"

    # Calculate timing
    SAMPLE_END=$(date +%s)
    SAMPLE_DURATION=$((SAMPLE_END - SAMPLE_START))
    PROCESSED=$((PROCESSED + 1))

    # Update running average
    if [ $PROCESSED -eq 1 ]; then
        AVG_TIME=$SAMPLE_DURATION
    else
        AVG_TIME=$(( (AVG_TIME * (PROCESSED - 1) + SAMPLE_DURATION) / PROCESSED ))
    fi

    # Calculate ETA
    REMAINING=$((TOTAL - COUNT))
    ETA_SECS=$((AVG_TIME * REMAINING))
    ETA_MINS=$((ETA_SECS / 60))

    # Update progress bar with completion
    progress_bar $COUNT $TOTAL "$SAMPLE" "✓ ${SAMPLE_DURATION}s" "ETA: ${ETA_MINS}m"
    echo ""
done

END_TIME=$(date +%s)
TOTAL_DURATION=$((END_TIME - START_TIME))
TOTAL_MINS=$((TOTAL_DURATION / 60))
TOTAL_SECS=$((TOTAL_DURATION % 60))

echo ""
echo "=============================================="
echo "All samples completed!"
echo "=============================================="
echo "Total time: ${TOTAL_MINS}m ${TOTAL_SECS}s"
echo "Processed:  $PROCESSED samples"
echo "Skipped:    $SKIPPED samples (already done)"
if [ $PROCESSED -gt 0 ]; then
    echo "Avg/sample: ${AVG_TIME}s"
fi
echo ""
echo "Output locations:"
echo "  Non-rRNA: 04_rRNA_filtered/non_rRNA/"
echo "  rRNA:     04_rRNA_filtered/rRNA/"
echo "  Logs:     logs/*_sortmerna.log"
