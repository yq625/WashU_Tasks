#!/bin/bash

REF="human_g1k_v37_decoy.fasta"

# Set the directories for each sample
SAMPLES=("Sample1" "Sample2" "Sample3")

# Run HaplotypeCaller for each sample
for SAMPLE in "${SAMPLES[@]}"; do
    gatk HaplotypeCaller \
        -R $REF \
        -I "${SAMPLE}/${SAMPLE}.aligned.sorted.bam" \
        -O "${SAMPLE}/${SAMPLE}.g.vcf.gz" \
        -ERC GVCF
done

# Combine GVCFs
COMBINED_GVCF="combined.g.vcf.gz"
gatk CombineGVCFs \
    -R $REF \
    --variant Sample1/Sample1.g.vcf.gz \
    --variant Sample2/Sample2.g.vcf.gz \
    --variant Sample3/Sample3.g.vcf.gz \
    -O $COMBINED_GVCF

# Genotype GVCFs
JOINT_VCF="joint.vcf.gz"
gatk GenotypeGVCFs \
    -R $REF \
    --variant $COMBINED_GVCF \
    -O $JOINT_VCF

# Apply VQSR for SNPs
RECAL_SNP="output.recal"
TRANCHES_SNP="output.tranches"
gatk VariantRecalibrator \
    -R $REF \
    -V $JOINT_VCF \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.vcf.gz \
    --resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.b37.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.b37.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_138.b37.vcf.gz \
    -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum \
    -mode SNP \
    -O $RECAL_SNP \
    --tranches-file $TRANCHES_SNP

gatk ApplyVQSR \
    -R $REF \
    -V $JOINT_VCF \
    --ts_filter_level 99.0 \
    --recal-file $RECAL_SNP \
    --tranches-file $TRANCHES_SNP \
    -mode SNP \
    -O "joint.filtered.vcf.gz"

echo "Task 07 pipeline complete."



