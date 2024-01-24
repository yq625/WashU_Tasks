#!/bin/bash

# Path to the reference genome
REF_GENOME="hg19.fa"

# Loop through each sample directory
for SAMPLE_DIR in Sample1 Sample2 Sample3
do
    echo "Processing $SAMPLE_DIR..."

    # Define file paths
    R1="$SAMPLE_DIR/${SAMPLE_DIR}.r1.fq.gz"
    R2="$SAMPLE_DIR/${SAMPLE_DIR}.r2.fq.gz"
    OUTPUT_SAM="$SAMPLE_DIR/${SAMPLE_DIR}.aligned.sam"

    # Align the reads to the reference genome
    bwa mem $REF_GENOME $R1 $R2 > $OUTPUT_SAM

    echo "$SAMPLE_DIR alignment complete."
done

echo "All samples processed."


