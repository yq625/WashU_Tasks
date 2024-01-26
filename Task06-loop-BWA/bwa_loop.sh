#!/bin/bash

# Path to the reference genome
REF_GENOME="human_g1k_v37_decoy.fasta"

# Loop through each sample directory
for SAMPLE_DIR in Sample1 Sample2 Sample3
do
    echo "Processing $SAMPLE_DIR..."

    # Define file paths
    R1="$SAMPLE_DIR/${SAMPLE_DIR}.r1.fq.gz"
    R2="$SAMPLE_DIR/${SAMPLE_DIR}.r2.fq.gz"
    RGFILE="$SAMPLE_DIR/${SAMPLE_DIR}.rgfile"

    OUTPUT_SAM="$SAMPLE_DIR/${SAMPLE_DIR}.aligned.sam"

    # Read the read group information from the rgfile
    RG=$(cat "$RGFILE")
    
    # Align the reads to the reference genome and add read group information
    bwa mem -M -R "$RG" "$REF_GENOME" "$R1" "$R2" > "$OUTPUT_SAM"


    echo "$SAMPLE_DIR alignment complete."
done

echo "All samples processed."


