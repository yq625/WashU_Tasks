# Task6: Whole Exome Sequencing Alignment Workflow

## Overview

This repository documents the workflow for aligning Whole Exome Sequencing data from three samples (Sample1, Sample2, and Sample3) to the human reference genome (GRCh37/hg19). The process involves downloading the reference genome, indexing it, and then aligning the sequencing data using the BWA (Burrows-Wheeler Aligner) tool.

## Downloading the Reference Genome

The GRCh37 reference genome was downloaded using `wget`. This version of the human genome is widely used for genomic studies.

### Steps:

1. **Download the Reference Genome**:
   ```bash
   wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
   ```
2. **Decompress the Genome File**:
   ```bash
   gunzip hg19.fa.gz
   ```

## Indexing the Reference Genome
Before alignment, the reference genome needs to be indexed using BWA.
### Command:
```bash
bwa index hg19.fa
```

## Aligning the Sequencing Data

The alignment of sequencing data for each sample (Sample1, Sample2, and Sample3) was conducted using the BWA MEM algorithm. The process is automated through a bash script named `bwa_loop.sh`.

### Execution:

To run the alignment script:
```bash
chmod +x bwa_loop.sh
./bwa_loop.sh
```
This script processes each sample directory, aligning the sequencing data to the GRCh37/hg19 reference genome and generating SAM files as output.

## Post-Alignment Processing
After alignment, the SAM files were further processed to convert them to sorted and indexed BAM files, which are more efficient for downstream analysis. The post-alignment processing is automated through `post_alignment_processing.sh`.

### Execution:
```bash
chmod +x post_alignment_processing.sh
./post_alignment_processing.sh
```
This script will iterate through each sample directory, convert the SAM file to BAM, sort the BAM file, and then index the sorted BAM file. It assumes that the SAM files are named in the format ${SAMPLE_DIR}.aligned.sam and are located in their respective sample directories. Make sure to replace Sample1, Sample2, and Sample3 with the actual names of your sample directories if they are different.

