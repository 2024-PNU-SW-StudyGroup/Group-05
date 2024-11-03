#!/bin/bash

# NGS Mapping Pipeline using BWA-MEM2
# Usage: ./mapping_pipeline.sh <reference> <read1> <read2> <output_prefix> <threads>

set -e

# Function to check command existence
check_command() {
    if ! command -v $1 &> /dev/null; then
        echo "Error: $1 is not installed"
        exit 1
    fi
}

# Check arguments
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <reference> <read1> <read2> <output_prefix> <threads>"
    exit 1
fi

# Assign arguments to variables
REFERENCE=$1
READ1=$2
READ2=$3
PREFIX=$4
THREADS=$5

# Check required tools
check_command bwa-mem2
check_command samtools
check_command java

# Create output directory
mkdir -p "${PREFIX}_output"
cd "${PREFIX}_output"

# Create log file
LOG="${PREFIX}.log"
echo "Starting mapping pipeline at $(date)" > "$LOG"

# Function to log commands and their outputs
run_cmd() {
    echo "Running: $@" >> "$LOG"
    "$@" 2>> "$LOG"
}

# Index reference if needed
if [ ! -f "${REFERENCE}.bwt.2bit.64" ]; then
    echo "Indexing reference genome..." | tee -a "$LOG"
    run_cmd bwa-mem2 index "$REFERENCE"
fi

# Mapping
echo "Performing mapping..." | tee -a "$LOG"
run_cmd bwa-mem2 mem -t "$THREADS" \
    -R "@RG\tID:${PREFIX}\tSM:${PREFIX}\tPL:ILLUMINA" \
    "$REFERENCE" "$READ1" "$READ2" \
    | samtools view -bh - > "${PREFIX}.bam"

# Sort BAM
echo "Sorting BAM..." | tee -a "$LOG"
run_cmd samtools sort -@ "$THREADS" \
    -o "${PREFIX}.sorted.bam" \
    "${PREFIX}.bam"

# Mark duplicates
echo "Marking duplicates..." | tee -a "$LOG"
run_cmd java -jar picard.jar MarkDuplicates \
    I="${PREFIX}.sorted.bam" \
    O="${PREFIX}.sorted.marked.bam" \
    M="${PREFIX}.marked_dup_metrics.txt"

# Index final BAM
echo "Indexing final BAM..." | tee -a "$LOG"
run_cmd samtools index "${PREFIX}.sorted.marked.bam"

# Generate mapping statistics
echo "Generating statistics..." | tee -a "$LOG"

# Basic statistics
run_cmd samtools flagstat "${PREFIX}.sorted.marked.bam" \
    > "${PREFIX}.flagstat.txt"

# Depth statistics
run_cmd samtools depth "${PREFIX}.sorted.marked.bam" \
    | awk '{sum+=$3} END {print "Average depth: "sum/NR}' \
    > "${PREFIX}.depth.txt"

# Insert size statistics
run_cmd samtools stats "${PREFIX}.sorted.marked.bam" \
    | grep "^IS" \
    > "${PREFIX}.insert_size.txt"

# Create summary report
cat << EOF > "${PREFIX}_summary.txt"
Mapping Summary Report
=====================
Date: $(date)
Sample: ${PREFIX}

Input Files:
- Reference: ${REFERENCE}
- Read1: ${READ1}
- Read2: ${READ2}

Mapping Statistics:
$(cat "${PREFIX}.flagstat.txt")

Coverage Statistics:
$(cat "${PREFIX}.depth.txt")

Duplicate Rate:
$(grep "PERCENT_DUPLICATION" "${PREFIX}.marked_dup_metrics.txt" | tail -n 1)

For detailed information, check:
- Duplicate metrics: ${PREFIX}.marked_dup_metrics.txt
- Insert size distribution: ${PREFIX}.insert_size.txt
- Full mapping statistics: ${PREFIX}.flagstat.txt
EOF

# Cleanup
rm "${PREFIX}.bam"

echo "Pipeline completed successfully!" | tee -a "$LOG"
echo "Check ${PREFIX}_summary.txt for results"