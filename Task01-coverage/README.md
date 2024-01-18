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
- `task1.chr22.bed`: BED file specifying regions of interest.

The results of the analysis will be stored in the `output` subdirectory within `Task1`.

## Usage Instructions
1. Navigate to the `Task1` directory.
2. Ensure that your BAM file (`task1.bam`) and BED file (`task1.chr22.bed`) are placed in this directory.
3. Run the scripts using the following commands:
   ```bash
   ./calculate_mean_coverage.sh
   ./calculate_coverage_10x_20x.sh
