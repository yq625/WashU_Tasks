# PLINK Data Merge: CHR22_AD and CHR22_CSF-Extension

## Project Overview

This project involves the merging of two genetic datasets, `CHR22_AD` and `CHR22_CSF-Extension`, which contain information on different individuals. The task aims to combine these datasets into a single set of files, ensuring proper alignment and integration of shared variants.

## Objective

The primary goal was to merge the `CHR22_AD` and `CHR22_CSF-Extension` datasets while addressing any variant mismatches. This ensures a comprehensive and accurate representation of the genetic information from both datasets.

## Methodology

### Initial Merge Attempt

The first step involved an attempt to merge the two datasets using PLINK. The command used was:

```bash
plink --bfile CHR22_AD --bmerge CHR22_CSF-Extension.bed CHR22_CSF-Extension.bim CHR22_CSF-Extension.fam --make-bed --out merged_attempt
```

### Handling Mismatches

The initial merge attempt revealed mismatches in variants, identified in the `merged_attempt-merge.missnp` file. To address this, the following step was taken:

**Exclude Mismatched SNPs**: Variants listed in the `.missnp` file were excluded from both datasets. The commands used were:

   ```bash
   plink --bfile CHR22_AD --exclude merged_attempt-merge.missnp --make-bed --out CHR22_AD_pruned
   plink --bfile CHR22_CSF-Extension --exclude merged_attempt-merge.missnp --make-bed --out CHR22_CSF-Extension_pruned
   ```

### Final Merge

After pruning the mismatched variants, the datasets were successfully merged using the following command:
```bash
plink --bfile CHR22_AD_pruned --bmerge CHR22_CSF-Extension_pruned.bed CHR22_CSF-Extension_pruned.bim CHR22_CSF-Extension_pruned.fam --make-bed --out final_merged_dataset
```
