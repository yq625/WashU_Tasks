#!/bin/bash

# Path to reference genome
REF_GENOME="path/to/hg19.fa"

# Loop through each sample directory
for SAMPLE_DIR in Sample1 Sample2 Sample3
do
    echo "Variant calling for $SAMPLE_DIR..."

    # Input recalibrated BAM and output GVCF file
    RECAL_BAM="$SAMPLE_DIR/${SAMPLE_DIR}.recal.bam"
    GVCF="$SAMPLE_DIR/${SAMPLE_DIR}.g.vcf.gz"

    # Call variants
    gatk HaplotypeCaller -R $REF_GENOME -I $RECAL_BAM -O $GVCF -ERC GVCF
done

