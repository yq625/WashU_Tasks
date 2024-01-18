#!/bin/bash

# Script to calculate mean coverage of regions in a BED file from a BAM file

# Ensure bedtools is installed
command -v bedtools >/dev/null 2>&1 || { echo >&2 "bedtools not found. Aborting."; exit 1; }

# Variables
BAM_FILE="task1.bam"
BED_FILE="task1.chr22.bed"
OUTPUT_DIR="output"
COVERAGE_FILE="${OUTPUT_DIR}/coverage.bedgraph"

# Create output directory if it doesn't exist
if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
fi

# Modify BED file to match the naming convention
sed 's/^chr//' $BED_FILE > "${OUTPUT_DIR}/modified_task1.chr22.bed"
BED_FILE="${OUTPUT_DIR}/modified_task1.chr22.bed"

# Generate a BEDGraph of coverage from the BAM file
bedtools genomecov -ibam $BAM_FILE -bg > $COVERAGE_FILE

# Intersect coverage with BED regions and calculate mean coverage
bedtools intersect -a $COVERAGE_FILE -b $BED_FILE | \
awk '{covSum += ($3 - $2) * $4; bases += ($3 - $2)} END {if (bases > 0) print "Mean Coverage: " covSum / bases; else print "No coverage data."}' > "${OUTPUT_DIR}/mean_coverage.txt"

