#!/bin/bash

# BLAST Analysis Pipeline
# Usage: ./blast_analysis.sh <query_fasta> <db_fasta> <output_dir>

# Exit on error
set -e

# Check arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <query_fasta> <db_fasta> <output_dir>"
    exit 1
fi

# Assign arguments to variables
QUERY=$1
DB_FASTA=$2
OUTPUT_DIR=$3

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check if BLAST+ is installed
if ! command_exists blastn; then
    echo "Error: BLAST+ is not installed"
    echo "Please install BLAST+ using:"
    echo "sudo apt-get install ncbi-blast+"
    echo "or"
    echo "conda install -c bioconda blast"
    exit 1
fi

# Create BLAST database
echo "Creating BLAST database..."
makeblastdb \
    -in "${DB_FASTA}" \
    -dbtype nucl \
    -out "${OUTPUT_DIR}/blastdb" \
    -parse_seqids

# Run BLAST with different output formats
echo "Running BLAST analysis..."

# Standard output
blastn \
    -query "${QUERY}" \
    -db "${OUTPUT_DIR}/blastdb" \
    -out "${OUTPUT_DIR}/blast_results.txt" \
    -evalue 1e-10 \
    -num_threads 4

# Tabular output
blastn \
    -query "${QUERY}" \
    -db "${OUTPUT_DIR}/blastdb" \
    -out "${OUTPUT_DIR}/blast_results.tsv" \
    -evalue 1e-10 \
    -num_threads 4 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

# XML output
blastn \
    -query "${QUERY}" \
    -db "${OUTPUT_DIR}/blastdb" \
    -out "${OUTPUT_DIR}/blast_results.xml" \
    -evalue 1e-10 \
    -num_threads 4 \
    -outfmt 5

# Create summary report
echo "Creating summary report..."
cat << EOF > "${OUTPUT_DIR}/analysis_report.txt"
BLAST Analysis Report
====================
Date: $(date)
Query: ${QUERY}
Database: ${DB_FASTA}
# 이전 스크립트에 이어서...

Parameters:
- E-value threshold: 1e-10
- Number of threads: 4
- Output formats: text, tabular, XML

Results Files:
- Standard output: ${OUTPUT_DIR}/blast_results.txt
- Tabular output: ${OUTPUT_DIR}/blast_results.tsv
- XML output: ${OUTPUT_DIR}/blast_results.xml

Summary Statistics:
$(grep -c "^>" "${OUTPUT_DIR}/blast_results.txt") total hits found
EOF

# Parse tabular results for basic statistics
if [ -f "${OUTPUT_DIR}/blast_results.tsv" ]; then
    echo -e "\nAlignment Statistics:" >> "${OUTPUT_DIR}/analysis_report.txt"
    echo "Average percent identity: $(awk '{sum+=$3} END {print sum/NR}' "${OUTPUT_DIR}/blast_results.tsv")" >> "${OUTPUT_DIR}/analysis_report.txt"
    echo "Average alignment length: $(awk '{sum+=$4} END {print sum/NR}' "${OUTPUT_DIR}/blast_results.tsv")" >> "${OUTPUT_DIR}/analysis_report.txt"
    echo "Average E-value: $(awk '{sum+=$11} END {print sum/NR}' "${OUTPUT_DIR}/blast_results.tsv")" >> "${OUTPUT_DIR}/analysis_report.txt"
fi

# Create visualization script
cat << 'EOF' > "${OUTPUT_DIR}/plot_results.R"
#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)

# Read BLAST results
data <- read.table("blast_results.tsv", 
                  sep="\t", 
                  col.names=c("qseqid","sseqid","pident","length",
                            "mismatch","gapopen","qstart","qend",
                            "sstart","send","evalue","bitscore"))

# Create identity distribution plot
pdf("identity_distribution.pdf")
ggplot(data, aes(x=pident)) +
    geom_histogram(bins=50, fill="steelblue", color="black") +
    theme_minimal() +
    labs(title="Distribution of Sequence Identity",
         x="Percent Identity",
         y="Count")
dev.off()

# Create alignment length distribution plot
pdf("length_distribution.pdf")
ggplot(data, aes(x=length)) +
    geom_histogram(bins=50, fill="darkgreen", color="black") +
    theme_minimal() +
    labs(title="Distribution of Alignment Lengths",
         x="Alignment Length",
         y="Count")
dev.off()

# Create e-value distribution plot
pdf("evalue_distribution.pdf")
ggplot(data, aes(x=-log10(evalue))) +
    geom_histogram(bins=50, fill="darkred", color="black") +
    theme_minimal() +
    labs(title="Distribution of E-values",
         x="-log10(E-value)",
         y="Count")
dev.off()
EOF

# Make R script executable
chmod +x "${OUTPUT_DIR}/plot_results.R"

# Run R script if R is installed
if command_exists Rscript; then
    echo "Generating plots..."
    cd "${OUTPUT_DIR}" && ./plot_results.R
    cd - > /dev/null
else
    echo "R is not installed. Skipping plot generation."
fi

echo "Analysis complete! Check ${OUTPUT_DIR} for results."