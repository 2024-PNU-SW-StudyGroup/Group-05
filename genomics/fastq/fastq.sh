#!/bin/bash

# FastQ Quality Control and Preprocessing Pipeline
# Usage: ./fastq_preprocess.sh <input_directory> <output_directory> <adapter_sequence>

# Exit on error
set -e

# Check arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_directory> <output_directory> <adapter_sequence>"
    echo "Example: $0 raw_data/ processed_data/ AGATCGGAAGAG"
    exit 1
fi

# Assign arguments to variables
INPUT_DIR=$1
OUTPUT_DIR=$2
ADAPTER=$3

# Create output directory structure
mkdir -p "${OUTPUT_DIR}/fastqc_raw"
mkdir -p "${OUTPUT_DIR}/fastqc_processed"
mkdir -p "${OUTPUT_DIR}/trimmed"

# Function to run FastQC
run_fastqc() {
    local input_dir=$1
    local output_dir=$2
    echo "Running FastQC on files in ${input_dir}"
    fastqc -t 4 -o "${output_dir}" "${input_dir}"/*.fastq*
}

# Function to run Cutadapt
run_cutadapt() {
    local input_file=$1
    local output_dir=$2
    local adapter=$3
    
    # Get filename without path and extension
    local basename=$(basename "${input_file}" .fastq)
    basename=$(basename "${basename}" .fq)
    basename=$(basename "${basename}" .gz)
    
    echo "Processing ${basename} with Cutadapt"
    
    # Run Cutadapt with quality trimming and adapter removal
    cutadapt \
        -a "${adapter}" \
        -q 20 \
        -m 50 \
        -o "${output_dir}/trimmed/${basename}_trimmed.fastq" \
        "${input_file}" \
        > "${output_dir}/trimmed/${basename}_cutadapt_report.txt"
}

# Main pipeline
echo "Starting QC and preprocessing pipeline..."
echo "Input directory: ${INPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Adapter sequence: ${ADAPTER}"

# Step 1: Initial FastQC
echo "Step 1: Running initial FastQC..."
run_fastqc "${INPUT_DIR}" "${OUTPUT_DIR}/fastqc_raw"

# Step 2: Cutadapt processing
echo "Step 2: Running Cutadapt..."
for file in "${INPUT_DIR}"/*.fastq*; do
    run_cutadapt "${file}" "${OUTPUT_DIR}" "${ADAPTER}"
done

# Step 3: FastQC on processed files
echo "Step 3: Running FastQC on processed files..."
run_fastqc "${OUTPUT_DIR}/trimmed" "${OUTPUT_DIR}/fastqc_processed"

# Create processing report
echo "Creating processing report..."
cat << EOF > "${OUTPUT_DIR}/processing_report.txt"
FastQ Processing Report
=====================
Date: $(date)
Input Directory: ${INPUT_DIR}
Output Directory: ${OUTPUT_DIR}
Adapter Sequence: ${ADAPTER}

Processed Files:
$(ls -1 "${INPUT_DIR}"/*.fastq*)

Quality Control Reports:
Raw Data FastQC: ${OUTPUT_DIR}/fastqc_raw
Processed Data FastQC: ${OUTPUT_DIR}/fastqc_processed

Cutadapt Reports: ${OUTPUT_DIR}/trimmed/*_cutadapt_report.txt
EOF

echo "Pipeline completed successfully!"
echo "Check ${OUTPUT_DIR}/processing_report.txt for details"