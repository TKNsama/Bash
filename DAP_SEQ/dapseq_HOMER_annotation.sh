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

module load HOMER
# Directory setup
ref_genome="/filer-5/agruppen/PBP/tan/indexDir/barley_index/220830_Bowman_pseudomolecules_and_unplaced_contigs_CPclean.fasta"
gtf="/filer-5/agruppen/PBP/tan/indexDir/barley_index/Bowman.gtf"
input_dir="/filer-5/user/tan/transfer/DAPseq/peaks"
output_dir="/filer-5/user/tan/transfer/DAPseq/peaks_analysis"

# 创建输出文件夹
mkdir -p "$output_dir"

# 循环遍历所有 .narrowPeak 文件
for peak_file in "$input_dir"/*.narrowPeak; do
    # 获取文件的基本名称，不带路径和扩展名
    base_name=$(basename "$peak_file" .narrowPeak)

    # 创建一个对应的输出文件夹
    mkdir -p "$output_dir/$base_name"

    # 1. Peak 注释
    echo "Annotating peaks for $base_name"
    annotatePeaks.pl "$peak_file" $ref_genome -gtf $gtf > "$output_dir/$base_name/annotated_peaks.txt"

    # 2. Motif 分析
    echo "Performing motif analysis for $base_name"
    findMotifsGenome.pl "$peak_file" $ref_genome "$output_dir/$base_name/motifs_output"

done

echo "HOMER analysis completed for all peaks."
