# Task 07: GATK Variant Calling Pipeline

## Overview
Task 07 involves processing genomic data through the GATK pipeline to perform variant calling. This task uses the hg19/b37 human reference genome and involves steps from initial variant calling to Variant Quality Score Recalibration (VQSR).

## Files and Directories
- `Sample1`, `Sample2`, `Sample3`: Directories containing aligned and sorted BAM files for each sample.
- VCF files for VQSR (`dbsnp_138.b37.vcf.gz`, `hapmap_3.3.b37.vcf.gz`, etc.).
- Reference genome files (`human_g1k_v37_decoy.fasta` and related index files).

## Pipeline Steps
1. **HaplotypeCaller**: Generates GVCFs for each sample.
2. **CombineGVCFs**: Combines individual GVCFs into a single GVCF.
3. **GenotypeGVCFs**: Performs joint genotyping to produce a raw VCF file.
4. **VariantRecalibrator and ApplyVQSR**: Applies VQSR using multiple resources like dbSNP, HapMap, and the 1000 Genomes Project.

## Step-by-Step Pipeline

### Step 1: Prepare the Reference Genome and VCF Files
Save reference genome and related index files given by the Task. Unzip them. Download vcf files:

```bash
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.indels.b37.vcf.gz
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/hapmap_3.3.b37.vcf.gz
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_omni2.5.b37.vcf.gz
```

### Step 2: Run Pipeline with `call_variant.sh`:

```bash
chmod +x call_variant.sh
./call_variant.sh
```


## Final Output

The pipeline's end product is a joint VCF file (`joint.vcf.gz`), which contains the consolidated variant calls from all three samples. This file is suitable for downstream genomic analyses and interpretations.

## Additional Information
- The script assumes the presence of necessary VCF files for VQSR in the working directory.
- For detailed information on GATK and variant calling, visit the [GATK website](https://gatk.broadinstitute.org/).
  
