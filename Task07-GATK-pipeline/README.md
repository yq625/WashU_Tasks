# Task 07: Variant Calling Pipeline for Whole Exome Sequencing

## Overview
This README outlines the specific steps for generating a joint VCF file from the Whole Exome Sequencing data of three samples. The process includes downloading the necessary reference files, processing BAM files, and performing variant calling using GATK v4.
### Required Input and Folder Structure
- **Reference Genome**: A FASTA file of the GRCh37/hg19 reference genome. 
- **dbSNP File**: A VCF file containing dbSNP variant annotations compatible with GRCh37/hg19.
- **Sample BAM Files**: Sorted and indexed BAM files for each of the three samples, ideally located in their respective directories (Sample1, Sample2, Sample3).
- **Scripts**: Ensure `process_bam.sh` and `variant_calling.sh` scripts are present in the working directory.
- **Conda Environment**: A Conda environment with GATK v4 installed and activated.

## Step-by-Step Pipeline

### Step 1: Download dbSNP File
```bash
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz
```
This file contains known variant annotations, which are essential for base quality score recalibration.

### Step 2: Process BAM Files
A script named `process_bam.sh` is used to mark duplicates and perform base quality score recalibration on each sample's BAM file.

To execute the script:
```bash
chmod +x process_bam.sh
./process_bam.sh
```

### Step 3: Variant Calling
Use the `variant_calling.sh` script to call variants on each sample's recalibrated BAM file using GATK's HaplotypeCaller.

To execute the script:
```bash
chmod +x variant_calling.sh
./variant_calling.sh
```

### Step 4: Joint Genotyping
After generating GVCFs for each sample, perform joint genotyping:
```bash
gatk CombineGVCFs -R $REF_GENOME --variant Sample1/Sample1.g.vcf.gz --variant Sample2/Sample2.g.vcf.gz --variant Sample3/Sample3.g.vcf.gz -O combined.g.vcf.gz
gatk GenotypeGVCFs -R $REF_GENOME --variant combined.g.vcf.gz -O output.vcf.gz
```

### Step 5: Post-Processing of VCF
Filter and refine the VCF file:
```bash
gatk VariantFiltration -R $REF_GENOME -V output.vcf.gz -O final_output.vcf.gz
```

## Final Output

The pipeline's end product is a joint VCF file (`final_joint_output.vcf.gz`), which contains the consolidated variant calls from all three samples. This file is suitable for downstream genomic analyses and interpretations.

## Conclusion

Following this pipeline, you will obtain a comprehensive view of the genetic variants present across the three exome-sequenced samples, facilitating further genomic studies and insights.
