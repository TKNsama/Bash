#!/bin/bash

#SBATCH --job-name=rnaSeqStep1Fastp                                 # Job name, will show up in squeue output
#SBATCH --auks=yes                                      # use Kerberos
#SBATCH --partition=cpu                                 # possible values: cpu, gpu
#SBATCH --ntasks=1                                      # Number of tasks (default=1)
#SBATCH --cpus-per-task=32                               # number of cpus for this task (default=1)
#SBATCH --mem=400G                                       # size of memory (default 10G in cpu and gpu)
#SBATCH --nodes=1                                       # Ensure that all cores are on one machine
#SBATCH --output=job_%j.out                             # File to which standard out will be written
#SBATCH --error=job_%j.err                              # File to which standard err will be written
#SBATCH --mail-type=END                                 # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=tan@ipk-gatersleben.de  # Email to which notifications will be sent
module load bwa
module load samtools
module load bcftools
module load seqkit

REF="Barley_Morex_V2_pseudomolecules.fa"
THREADS=32
K=1000  # 分片大小
INPUT_DIR="/filer-5/user/tan/pangenome/selected_library_files/Barley/Barley_Pangenome_v2_2023/Assemblies/compressed_assemblies"
RESULT_DIR="results"
mkdir -p "$RESULT_DIR"

# 创建 BWA 索引（只需一次）
if [[ ! -f "${REF}.bwt" ]]; then
  echo "Building BWA index for $REF..."
  bwa index "$REF"
fi

for fasta in "$INPUT_DIR"/*.fasta; do
  base=$(basename "$fasta" .fasta)
  echo "Processing $base..."

  SPLIT_FASTA="$RESULT_DIR/$base.split.fa"
  BAM="$RESULT_DIR/$base.bam"
  VCF="$RESULT_DIR/$base.snps.vcf"

  # 如果VCF已存在，跳过当前样本
  if [[ -f "$VCF" ]]; then
    echo "VCF file $VCF already exists. Skipping $base."
    continue
  fi

  # 1. 拆分 fasta 为固定 K 窗口
  if [[ ! -f "$SPLIT_FASTA" ]]; then
    echo "Splitting $fasta into $K bp windows..."
    seqkit sliding -s "$K" -W "$K" "$fasta" -o "$SPLIT_FASTA"
  fi

  # 2. 比对到参考基因组
  if [[ ! -f "$BAM" ]]; then
    echo "Aligning to reference with BWA..."
    bwa mem -t "$THREADS" "$REF" "$SPLIT_FASTA" | samtools sort -@ "$THREADS" -o "$BAM"
    samtools index "$BAM"
  fi

  # 3. SNP/Indel calling
  echo "Calling SNPs with bcftools..."
  bcftools mpileup -Ou -f "$REF" "$BAM" | bcftools call -mv -Ov -o "$VCF"

done

