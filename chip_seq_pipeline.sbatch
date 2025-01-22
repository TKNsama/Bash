#!/bin/bash

#SBATCH --job-name=rnaSeqStep1Fastp
#SBATCH --auks=yes
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=tan@ipk-gatersleben.de
#SBATCH --time=168:00:00

# Load necessary modules
module load fastp
module load bowtie2
module load samtools
module load macs
module load bedtools
module load bedGraphToBigWig

# 参数设置
RAW_DATA_DIR="/filer-5/agruppen/PBP/tan/chipseq"                # 原始数据目录
workspace="/filer-5/agruppen/PBP/tan/chipseq/workspace"                    # 输出目录
REFERENCE_GENOME="/filer-5/agruppen/PBP/tan/indexDir/barley_index/Barley_Morex_V2_pseudomolecules_bowtie2_index"       # Bowtie2 索引的基因组前缀
SAMPLE_NAME="SRR1223012_1"                     # 样本名称
INPUT_NAME="SRR1223527_1"                       # 对照样本名称
THREADS=8                                # 使用的线程数

# 创建输出目录
mkdir -p ${workspace}/fastp
mkdir -p ${workspace}/bowtie2
mkdir -p ${workspace}/macs2

echo "### Step 1: 数据质量控制与修剪 ###"
# 使用 fastp 进行质量控制和修剪
fastp \
  -i ${RAW_DATA_DIR}/${SAMPLE_NAME}.fastq.gz \
  -o ${workspace}/fastp/${SAMPLE_NAME}_clean.fastq.gz \
  -h ${workspace}/fastp/${SAMPLE_NAME}_fastp.html \
  -j ${workspace}/fastp/${SAMPLE_NAME}_fastp.json \
  --thread ${THREADS}

fastp \
  -i ${RAW_DATA_DIR}/${INPUT_NAME}.fastq.gz \
  -o ${workspace}/fastp/${INPUT_NAME}_clean.fastq.gz \
  -h ${workspace}/fastp/${INPUT_NAME}_fastp.html \
  -j ${workspace}/fastp/${INPUT_NAME}_fastp.json \
  --thread ${THREADS}

echo "### Step 2: 读取比对 ###"
# 使用 bowtie2 将读取比对到参考基因组
bowtie2 \
  -x ${REFERENCE_GENOME} \
  -U ${workspace}/fastp/${SAMPLE_NAME}_clean.fastq.gz \
  -S ${workspace}/bowtie2/${SAMPLE_NAME}.sam \
  -p ${THREADS}

bowtie2 \
  -x ${REFERENCE_GENOME} \
  -U ${workspace}/fastp/${INPUT_NAME}_clean.fastq.gz \
  -S ${workspace}/bowtie2/${INPUT_NAME}.sam \
  -p ${THREADS}

echo "### Step 3: SAM 文件转换为 BAM 文件并去重 ###"
# 将 SAM 文件转换为 BAM 文件并进行排序和去重
samtools view -bS ${workspace}/bowtie2/${SAMPLE_NAME}.sam | \
  samtools sort -o ${workspace}/bowtie2/${SAMPLE_NAME}_sorted.bam

samtools view -bS ${workspace}/bowtie2/${INPUT_NAME}.sam | \
  samtools sort -o ${workspace}/bowtie2/${INPUT_NAME}_sorted.bam

# 使用 samtools markdup 去除 PCR 冗余
samtools markdup -r ${workspace}/bowtie2/${SAMPLE_NAME}_sorted.bam ${workspace}/bowtie2/${SAMPLE_NAME}_dedup.bam
samtools markdup -r ${workspace}/bowtie2/${INPUT_NAME}_sorted.bam ${workspace}/bowtie2/${INPUT_NAME}_dedup.bam

# 索引 BAM 文件
samtools index -c ${workspace}/bowtie2/${SAMPLE_NAME}_dedup.bam
samtools index -c ${workspace}/bowtie2/${INPUT_NAME}_dedup.bam

echo "### Step 4: 峰值调用 ###"
# 使用 MACS2 进行峰值调用
macs3 callpeak \
  -t ${workspace}/bowtie2/${SAMPLE_NAME}_dedup.bam \
  -c ${workspace}/bowtie2/${INPUT_NAME}_dedup.bam \
  -f BAM \
  -g hs \
  -n ${SAMPLE_NAME} \
  --outdir ${workspace}/macs2 \
  --nomodel \
  --extsize 200 \
  -q 0.05

samtools depth ${workspace}/bowtie2/${SAMPLE_NAME}_dedup.bam > ${SAMPLE_NAME}_depth.txt
samtools depth ${workspace}/bowtie2/${INPUT_NAME}_dedup.bam > ${INPUT_NAME}_depth.txt

bedtools genomecov -ibam ${workspace}/bowtie2/${SAMPLE_NAME}_dedup.bam -bg > ${SAMPLE_NAME}_depth.bedgraph
bedtools genomecov -ibam ${workspace}/bowtie2/${INPUT_NAME}_dedup.bam -bg > ${INPUT_NAME}_depth.bedgraph

bedGraphToBigWig ${SAMPLE_NAME}_depth.bedgraph /filer-5/agruppen/PBP/tan/indexDir/barley_index/Barley_Morex_V2_pseudomolecules.chrom.sizes ${SAMPLE_NAME}.bw
bedGraphToBigWig ${INPUT_NAME}_depth.bedgraph /filer-5/agruppen/PBP/tan/indexDir/barley_index/Barley_Morex_V2_pseudomolecules.chrom.sizes ${INPUT_NAME}.bw


echo "### Pipeline 完成 ###"



