#!/bin/bash

# ibd_analysis.sh

# Function to check if PLINK is installed
function check_plink_installed() {
    if ! command -v plink &> /dev/null
    then
        echo "PLINK could not be found. Please install PLINK to continue."
        exit 1
    else
        echo "PLINK is installed."
    fi
}

# Function to convert VCF to PLINK format and perform IBD analysis
function perform_ibd_analysis() {
    local vcf_file="task3.vcf.gz"
    
    # Check if the input VCF file exists
    if [ ! -f "$vcf_file" ]; then
        echo "VCF file does not exist: $vcf_file"
        exit 1
    fi

    # Decompress the VCF file if it is gzipped
    echo "Decompressing $vcf_file..."
    gunzip -k "$vcf_file"
    vcf_file="${vcf_file%.gz}"  # Update the vcf_file variable to the decompressed file

    # Convert VCF to PLINK binary format
    echo "Converting VCF to PLINK binary format..."
    plink --vcf "$vcf_file" --make-bed --out task3

    # Perform IBD analysis
    echo "Performing IBD analysis..."
    plink --bfile task3 --genome --out task3_ibd

    echo "IBD analysis complete. Output files are task3_ibd.genome and related binary files."
}

# Check if PLINK is installed
check_plink_installed

# Perform IBD analysis on the provided VCF file
perform_ibd_analysis

