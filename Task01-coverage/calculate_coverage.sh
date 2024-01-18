#!/bin/bash

# Ensure bedtools is installed
command -v bedtools >/dev/null 2>&1 || { echo >&2 "bedtools not found. Aborting."; exit 1; }

# Variables
BAM_FILE="task1.bam"
BED_FILE="task1.chr22.bed"
OUTPUT_DIR="output"
COVERAGE_FILE="${OUTPUT_DIR}/coverage.bedgraph"
MODIFIED_BED_FILE="${OUTPUT_DIR}/modified_task1.chr22.bed"

# Create output directory if it doesn't exist
if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
fi

# Modify BED file to match the naming convention
sed 's/^chr//' $BED_FILE > $MODIFIED_BED_FILE

# Generate a BEDGraph of coverage from the BAM file
bedtools genomecov -ibam $BAM_FILE -bg > $COVERAGE_FILE

# Intersect coverage with modified BED regions
bedtools intersect -a $COVERAGE_FILE -b $MODIFIED_BED_FILE > "${OUTPUT_DIR}/coverage_intersect.bedgraph"

# Calculate mean coverage
mean_coverage=$(awk '{covSum += ($3 - $2) * $4; bases += ($3 - $2)} END {if (bases > 0) print "Mean Coverage: " covSum / bases; else print "No coverage data."}' "${OUTPUT_DIR}/coverage_intersect.bedgraph")

# Calculate coverage at 10X and 20X
for COVERAGE in 10 20; do
    awk -v COV=$COVERAGE '$4 >= COV' "${OUTPUT_DIR}/coverage_intersect.bedgraph" | \
    awk '{total += ($3 - $2)} END {print total}' > "${OUTPUT_DIR}/coverage_${COVERAGE}x.txt"
done

# Calculate total length of regions
total_length=$(awk '{total += ($3 - $2)} END {print total}' $MODIFIED_BED_FILE)

# Calculate percentages
for COVERAGE in 10 20; do
    covered_length=$(cat "${OUTPUT_DIR}/coverage_${COVERAGE}x.txt")
    percentage=$(echo "scale=2; $covered_length / $total_length * 100" | bc)
    echo "Percentage Coverage at ${COVERAGE}X: $percentage" > "${OUTPUT_DIR}/percentage_${COVERAGE}x.txt"
done

# Create the final output file
echo "Mean Coverage: $mean_coverage" > "${OUTPUT_DIR}/coverage_summary.txt"
echo "Percentage Coverage at 10X: $(cat "${OUTPUT_DIR}/percentage_10x.txt")%" >> "${OUTPUT_DIR}/coverage_summary.txt"
echo "Percentage Coverage at 20X: $(cat "${OUTPUT_DIR}/percentage_20x.txt")%" >> "${OUTPUT_DIR}/coverage_summary.txt"

# Display results
cat "${OUTPUT_DIR}/coverage_summary.txt"

