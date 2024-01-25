#!/bin/bash

# Path to reference genome and dbSNP VCF
REF_GENOME="path/to/hg19.fa"
DBSNP="path/to/00-All.vcf.gz"

# Loop through each sample directory
for SAMPLE_DIR in Sample1 Sample2 Sample3
do
    echo "Processing $SAMPLE_DIR..."

    # Define input and output file paths
    INPUT_BAM="$SAMPLE_DIR/${SAMPLE_DIR}.bam"
    MARKED_BAM="$SAMPLE_DIR/${SAMPLE_DIR}.marked.bam"
    METRICS_FILE="$SAMPLE_DIR/${SAMPLE_DIR}.metrics.txt"
    RECAL_TABLE="$SAMPLE_DIR/${SAMPLE_DIR}.recal.table"
    RECAL_BAM="$SAMPLE_DIR/${SAMPLE_DIR}.recal.bam"

    # Mark duplicates
    gatk MarkDuplicates -I $INPUT_BAM -O $MARKED_BAM -M $METRICS_FILE

    # Base recalibration
    gatk BaseRecalibrator -R $REF_GENOME --known-sites $DBSNP -I $MARKED_BAM -O $RECAL_TABLE
    gatk ApplyBQSR -R $REF_GENOME -I $MARKED_BAM --bqsr-recal-file $RECAL_TABLE -O $RECAL_BAM
done

