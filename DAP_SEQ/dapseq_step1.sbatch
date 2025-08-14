#!/bin/bash

#SBATCH --job-name=DAPseq_Pipeline
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tan@ipk-gatersleben.de
# Load modules
module load fastp
module load bwa
module load samtools
module load MACS3

# Directory setup
raw_data_dir="/filer-5/user/tan/transfer/vrs5_dap_seq" # raw data location
clean_data_dir="/filer-5/user/tan/transfer/vrs5_dap_seq/clean"
bam_dir="/filer-5/user/tan/transfer/vrs5_dap_seq/bam"
peak_dir="/filer-5/user/tan/transfer/vrs5_dap_seq/peak"
ref_genome="/filer-5/agruppen/PBP/tan/indexDir/barley_index/Barley_Morex_V2_pseudomolecules.fa" # ref genome
threads=32

# Create output dirs
mkdir -p "$clean_data_dir" "$bam_dir" "$peak_dir"

# Loop over FASTQ files
for fq1 in "$raw_data_dir"/*_1.fastq.gz; do
    fq2="${fq1/_1.fastq.gz/_2.fastq.gz}"
    sample_name=$(basename "$fq1" _1.fastq.gz)

    echo "Processing sample: $sample_name"

    # Step 1: fastp for QC
    if [[ ! -f "$clean_data_dir/${sample_name}_clean_R1.fastq.gz" ]]; then
        echo "[FASTP] Running QC for $sample_name..."
        fastp -i "$fq1" -I "$fq2" \
              -o "$clean_data_dir/${sample_name}_clean_R1.fastq.gz" \
              -O "$clean_data_dir/${sample_name}_clean_R2.fastq.gz" \
              --thread $threads --detect_adapter_for_pe \
              --html "$clean_data_dir/${sample_name}_fastp.html"
    else
        echo "[FASTP] Skipping $sample_name (already exists)"
    fi

    # Step 2: Align with BWA
    if [[ ! -f "$bam_dir/${sample_name}.sorted.bam" ]]; then
        echo "[BWA] Aligning $sample_name..."
        bwa mem -t $threads "$ref_genome" \
            "$clean_data_dir/${sample_name}_clean_R1.fastq.gz" \
            "$clean_data_dir/${sample_name}_clean_R2.fastq.gz" \
            | samtools view -@ $threads -bS - \
            | samtools sort -@ $threads -o "$bam_dir/${sample_name}.sorted.bam"

        samtools index "$bam_dir/${sample_name}.sorted.bam"
    else
        echo "[BWA] Skipping $sample_name (BAM already exists)"
    fi

    # Step 3: Peak calling with MACS3
    if [[ ! -f "$peak_dir/${sample_name}_peaks.narrowPeak" ]]; then
        echo "[MACS3] Calling peaks for $sample_name..."
        macs3 callpeak -t "$bam_dir/${sample_name}.sorted.bam" \
            -f BAMPE -g 1.2e8 -n "$sample_name" \
            --outdir "$peak_dir" \
            --nomodel --shift -100 --extsize 200 \
            --keep-dup all -q 0.01
    else
        echo "[MACS3] Skipping $sample_name (peaks already exist)"
    fi

done
