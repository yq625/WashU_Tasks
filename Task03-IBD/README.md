# Task3 IBD Analysis

## Overview
This repository contains the scripts and results for the IBD (Identity By Descent) analysis for Task3. The project involves processing genetic data in `.vcf.gz` format and analyzing the relationships between individuals using PLINK.

## Contents
- `ibd_analysis.sh`: A shell script to perform IBD analysis using PLINK.
- `plot_pedigree.py`: A Python script to visualize the results of the IBD analysis.
- `pedigree_chart_full.png`: The output visual representation of the IBD analysis.
- `task3_ibd.genome`: The output file from PLINK with IBD analysis results.

## How to Run

### IBD Analysis
To perform the IBD analysis, make sure you have PLINK installed and the `task3.vcf.gz` file in the same directory as `ibd_analysis.sh`. Run the following command in your terminal:

```bash
./ibd_analysis.sh
```
This script checks for PLINK, converts the .vcf.gz file to PLINK format, and performs the IBD analysis.

### Visualization
Make sure Python is installed with the pandas, matplotlib, and networkx libraries.
Run the plot_pedigree.py script to generate a visualization:
```bash
python plot_pedigree.py
```
This script creates a plot of the relationships based on the PI_HAT values from the IBD analysis.

## Results
The pedigree_chart_full.png file in the repository shows the visual representation of the relationships identified in the IBD analysis. Lines in the plot indicate the strength of genetic relationships (PI_HAT value).

## Dependencies
PLINK (required for ibd_analysis.sh)
Python with libraries: pandas, matplotlib, networkx (required for plot_pedigree.py)
