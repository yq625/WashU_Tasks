# Task 1: BAM File Coverage Analysis

## Project Description
This repository contains scripts for analyzing the coverage depth of a BAM file in specific genomic regions defined in a BED file. The analysis includes calculating the mean coverage across these regions and determining the percentage of these regions covered at a depth of 10X and 20X.

## Prerequisites
To run these scripts, `bedtools` should be installed on your system. You can find installation instructions for `bedtools` [here](https://bedtools.readthedocs.io/en/latest/content/installation.html).

## Usage

1. Place your BAM file (e.g., `task1.bam`), BAI file (e.g., `task1.bai`), and BED file (e.g., `task1.chr22.bed`) in the same directory as this script (`calculate_coverage.sh`).

2. Execute the script using the following command:

   ```bash
   ./calculate_coverage.sh
## Output
Upon running the analysis scripts, the following output file will be generated in the `output` directory within the `Task1` folder:
- `coverage_summary.txt`: Contains the mean coverage and percentage coverage at 10X and 20X.
- `coverage.bedgraph`: BEDGraph file with coverage information.
- `modified_task1.chr22.bed`: Modified BED file.
- `coverage_intersect.bedgraph`: Intersection of coverage with modified BED regions.
- `coverage_10x.txt`: Coverage at 10X.
- `coverage_20x.txt`: Coverage at 20X.
- `percentage_10x.txt`: Percentage coverage at 10X.
- `percentage_20x.txt`: Percentage coverage at 20X.

You can access these files to review the calculated coverage and percentages.
