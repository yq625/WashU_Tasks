# WashU_Tasks

This repository contains scripts and documentation for various tasks in genomic data analysis, including coverage analysis, merging bfiles, and variant calling using BWA and GATK.

## Important Note on Data Files
Due to the large size of input and output files, especially for Tasks 4, 6, and 7, they are not available for download from this repository. The scripts and methodologies are provided, and users are encouraged to apply them to the datasets given by WashU.

## Repository Structure
- `Task01-coverage`: Scripts related to coverage analysis.
- `Task03-IBD`: Documentation and scripts for Identity by Descent analysis.
- `Task04-merge-bfiles`: Scripts for merging binary files in genomic datasets.
- `Task06-loop-BWA`: Script for aligning sequences using BWA in a loop across multiple samples.
- `Task07-GATK-pipeline`: Scripts for processing genomic data through the GATK pipeline, including variant calling and quality score recalibration.

## Workflow for Tasks 06 and 07
- Both Task 06 and Task 07 scripts operate in the same working directory.
- This directory contains `Sample1`, `Sample2`, and `Sample3` folders with necessary data files.
- Ensure that the working directory is properly set up with these folders before executing the scripts.

## Using the Scripts
- Each task folder contains a README with specific instructions for running the scripts.
- Ensure all prerequisites (like GATK, BWA) are installed and properly configured.
- Replace file paths in scripts with the correct locations on your system.

For more detailed information on each task and to access the scripts, please navigate to the respective directories in this repository.
