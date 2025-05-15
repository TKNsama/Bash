#!/bin/bash

#SBATCH --job-name=rnaSeqStep1Fastp                                 # Job name, will show up in squeue output
#SBATCH --auks=yes                                      # use Kerberos
#SBATCH --partition=cpu                                 # possible values: cpu, gpu
#SBATCH --ntasks=1                                      # Number of tasks (default=1)
#SBATCH --cpus-per-task=12                               # number of cpus for this task (default=1)
#SBATCH --mem=100G                                       # size of memory (default 10G in cpu and gpu)
#SBATCH --nodes=1                                       # Ensure that all cores are on one machine
#SBATCH --output=job_%j.out                             # File to which standard out will be written
#SBATCH --error=job_%j.err                              # File to which standard err will be written
#SBATCH --mail-type=END                                 # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=tan@ipk-gatersleben.de  # Email to which notifications will be sent
#SBATCH --time=18:16:40

# don't know why, but it should be inputed
source /etc/profile
module load fastp
module load bwa
module load samtools
module load macs
# Directory setup
raw_data_dir="/filer-5/user/tan/transfer/DAPseq"                  # <-- 修改为你的路径
clean_data_dir="/filer-5/user/tan/transfer/DAPseq/clean"
bam_dir="/filer-5/user/tan/transfer/DAPseq/bam"
peak_dir="/filer-5/user/tan/transfer/DAPseq/peaks"
ref_genome="/filer-5/agruppen/PBP/tan/indexDir/barley_index/220830_Bowman_pseudomolecules_and_unplaced_contigs_CPclean.fasta"
threads=12

# Create output dirs
mkdir -p "$clean_data_dir" "$bam_dir" "$peak_dir"

# Loop over FASTQ files
for fq in "$raw_data_dir"/*.fastq.gz; do
    sample_name=$(basename "$fq" .fastq.gz)

    echo "Processing sample: $sample_name"

    # Step 1: fastp for QC
    fastp -i "$fq" \
          -o "$clean_data_dir/${sample_name}_clean.fastq.gz" \
          --thread $threads --detect_adapter_for_pe \
          --html "$clean_data_dir/${sample_name}_fastp.html"

    # Step 2: Align with BWA
    bwa mem -t $threads "$ref_genome" \
        "$clean_data_dir/${sample_name}_clean.fastq.gz" \
        | samtools view -@ $threads -bS - \
        | samtools sort -@ $threads -o "$bam_dir/${sample_name}.sorted.bam"

    # Index BAM
    samtools index -c "$bam_dir/${sample_name}.sorted.bam"

    # Step 3: Peak calling with MACS2 (单端数据用 BAM 格式)
    macs callpeak -t "$bam_dir/${sample_name}.sorted.bam" \
        -f BAM -g 4.8e9 -n "$sample_name" \
        --outdir "$peak_dir" \
        --nomodel --shift -100 --extsize 200 \
        --keep-dup all -q 0.01
done

