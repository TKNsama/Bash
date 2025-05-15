#!/bin/bash

#SBATCH --job-name=DAPseq_Pipeline
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tan@ipk-gatersleben.de

# Load modules
source /etc/profile
module load fastp
module load bwa
module load samtools
module load macs2

# Directory setup
raw_data_dir="/path/to/raw_fastq"                  # <-- 修改为你的路径
clean_data_dir="/path/to/clean_fastq"
bam_dir="/path/to/bam"
peak_dir="/path/to/peaks"
ref_genome="/path/to/genome.fa"                    # <-- 修改为你的参考基因组路径
threads=12

# Create output dirs
mkdir -p "$clean_data_dir" "$bam_dir" "$peak_dir"

# Loop over FASTQ files
for fq1 in "$raw_data_dir"/*_R1.fastq.gz; do
    fq2="${fq1/_R1.fastq.gz/_R2.fastq.gz}"
    sample_name=$(basename "$fq1" _R1.fastq.gz)

    echo "Processing sample: $sample_name"

    # Step 1: fastp for QC
    fastp -i "$fq1" -I "$fq2" \
          -o "$clean_data_dir/${sample_name}_clean_R1.fastq.gz" \
          -O "$clean_data_dir/${sample_name}_clean_R2.fastq.gz" \
          --thread $threads --detect_adapter_for_pe \
          --html "$clean_data_dir/${sample_name}_fastp.html"

    # Step 2: Align with BWA
    bwa mem -t $threads "$ref_genome" \
        "$clean_data_dir/${sample_name}_clean_R1.fastq.gz" \
        "$clean_data_dir/${sample_name}_clean_R2.fastq.gz" \
        | samtools view -@ $threads -bS - \
        | samtools sort -@ $threads -o "$bam_dir/${sample_name}.sorted.bam"

    # Optional: mark duplicates or index BAM
    samtools index "$bam_dir/${sample_name}.sorted.bam"

    # Step 3: Peak calling with MACS2 (can adjust -f BAMPE or BAM as needed)
    macs2 callpeak -t "$bam_dir/${sample_name}.sorted.bam" \
        -f BAMPE -g 1.2e8 -n "$sample_name" \
        --outdir "$peak_dir" \
        --nomodel --shift -100 --extsize 200 \
        --keep-dup all -q 0.01
done
