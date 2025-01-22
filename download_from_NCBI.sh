#!/bin/bash

#SBATCH --job-name=rnaSeqStep1Fastp
#SBATCH --auks=yes
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=tan@ipk-gatersleben.de

# Load necessary modules
module load sratoolkit

# Define data directory and SRA list file
data_dir="/filer-5/agruppen/PBP/tan/chipseq"
sra_list="SRA_list.txt"

# Create data directory if it doesn't exist
mkdir -p "$data_dir"

# Iterate through each SRA accession in the SRA_list.txt file
while IFS= read -r sra_acc; do
    if [[ ! -z "$sra_acc" ]]; then
        echo "Processing $sra_acc..."
        # Download the SRA file
        prefetch --output-directory "$data_dir" "$sra_acc"
        
        # Define the path to the downloaded SRA file
        sra_file="$data_dir/$sra_acc.sra"
        
        # Check if the SRA file was successfully downloaded
        if [[ -f "$sra_file" ]]; then
            echo "Extracting $sra_file..."
            # Use fastq-dump to extract the SRA file to fastq format
            fastq-dump --gzip --split-files --outdir "$data_dir" "$sra_file"
        else
            echo "Failed to download $sra_acc"
        fi
    fi
done < "$sra_list"

echo "All tasks completed."

