#!/bin/bash

# Loop through each sample directory
for SAMPLE_DIR in Sample1 Sample2 Sample3
do
    echo "Processing post-alignment for $SAMPLE_DIR..."

    # Define input SAM and output BAM file paths
    INPUT_SAM="$SAMPLE_DIR/${SAMPLE_DIR}.aligned.sam"
    OUTPUT_BAM="$SAMPLE_DIR/${SAMPLE_DIR}.aligned.bam"
    SORTED_BAM="$SAMPLE_DIR/${SAMPLE_DIR}.aligned.sorted.bam"

    # Check if the input SAM file exists
    if [ ! -f $INPUT_SAM ]; then
        echo "Error: SAM file for $SAMPLE_DIR does not exist."
        continue
    fi

    # Convert SAM to BAM
    echo "Converting SAM to BAM for $SAMPLE_DIR..."
    samtools view -b $INPUT_SAM > $OUTPUT_BAM

    # Sort the BAM file
    echo "Sorting BAM file for $SAMPLE_DIR..."
    samtools sort $OUTPUT_BAM -o $SORTED_BAM

    # Index the sorted BAM file
    echo "Indexing the sorted BAM file for $SAMPLE_DIR..."
    samtools index $SORTED_BAM

    echo "$SAMPLE_DIR post-alignment processing complete."
done

echo "All samples processed."
