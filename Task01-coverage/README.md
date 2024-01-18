# Task 1: BAM File Coverage Analysis

## Project Description
This repository contains scripts for analyzing the coverage depth of a BAM file in specific genomic regions defined in a BED file. The analysis includes calculating the mean coverage across these regions and determining the percentage of these regions covered at a depth of 10X and 20X.

## Prerequisites
To run these scripts, `bedtools` should be installed on your system. You can find installation instructions for `bedtools` [here](https://bedtools.readthedocs.io/en/latest/content/installation.html).

## Structure of the Repository
The `Task1` directory contains both the scripts and the data files:
- `calculate_mean_coverage.sh`: Script to calculate the mean coverage over specified regions.
- `calculate_coverage_10x_20x.sh`: Script to calculate the percentage of regions covered at 10X and 20X.
- `task1.bam`: BAM file containing sequencing data.
- `task1.bai`: BAI file, the index of the BAM file.
- `task1.chr22.bed`: BED file specifying regions of interest.

The results of the analysis will be stored in the `output` subdirectory within `Task1`.

## Usage Instructions
1. Navigate to the `Task1` directory.
2. Ensure that your BAM file (`task1.bam`), BAI file (`task1.bai`), and BED file (`task1.chr22.bed`) are placed in this directory.
3. Make the shell scripts executable using the following commands:
   ```bash
   chmod +x calculate_mean_coverage.sh
   chmod +x calculate_coverage_10x_20x.sh
4. Run the scripts using the following commands:
   ```bash
   ./calculate_mean_coverage.sh
   ./calculate_coverage_10x_20x.sh
## Output
Upon running the analysis scripts, the following output files will be generated in the `output` directory within the `Task1` folder:

- `mean_coverage.txt`: This file contains the mean coverage across the specified genomic regions in the BED file. It provides an average coverage depth for the regions of interest.

- `percentage_10x.txt`: This file contains the percentage of genomic regions covered at a minimum depth of 10X. It helps assess the coverage at a specific threshold.

- `percentage_20x.txt`: Similar to the `percentage_10x.txt` file, this one contains the percentage of genomic regions covered at a minimum depth of 20X.

You can access and analyze these output files to gain insights into the coverage of your BAM file within the specified regions of interest.
