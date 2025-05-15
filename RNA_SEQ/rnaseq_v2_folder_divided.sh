#!/bin/bash

#SBATCH --job-name=rna_seq                                 # Job name, will show up in squeue output
#SBATCH --auks=yes                                      # use Kerberos
#SBATCH --partition=cpu                                 # possible values: cpu, gpu
#SBATCH --ntasks=1                                      # Number of tasks (default=1)
#SBATCH --cpus-per-task=8                               # number of cpus for this task (default=1)
#SBATCH --mem=100G                                       # size of memory (default 10G in cpu and gpu)
#SBATCH --nodes=1                                       # Ensure that all cores are on one machine
#SBATCH --output=job_%j.out                             # File to which standard out will be written
#SBATCH --error=job_%j.err                              # File to which standard err will be written
#SBATCH --mail-type=END                                 # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=tan@ipk-gatersleben.de  # Email to which notifications will be sent 

# rawdata需要后缀为_1.fq.gz _2.fq.gz

# ===== 模块加载 =====
module load fastp
module load kallisto

# ===== 参数设置 =====
raw_data_dir="/filer-5/agruppen/PBP/huang/MADS5/RNA_Seq/VRN1-HA/raw_data/X204SC23121144-Z01-F003/Merged"
processed_data_dir="/filer-5/agruppen/PBP/tan/rna_seq/mads5_2"
index_dir="/filer-5/agruppen/PBP/tan/indexDir"
index_file="$index_dir/barley_v2_ALL_cds.idx"
threads=8

# ===== 输出目录准备 =====
clean_dir="$processed_data_dir/kallisto"
mkdir -p "$clean_dir"

# ===== 日志函数 =====
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# ===== 主循环处理样本 =====
for sample_folder in "$raw_data_dir"/N*; do
    if [[ -d "$sample_folder" ]]; then
        fq1=$(ls "$sample_folder"/*_1.fq.gz 2>/dev/null | head -n 1)
        fq2=$(echo "$fq1" | sed 's/_1.fq.gz/_2.fq.gz/')
        
        if [[ -f "$fq1" && -f "$fq2" ]]; then
            sample_name=$(basename "$fq1" | sed 's/_1.fq.gz//')
            clean_fq1="${clean_dir}/${sample_name}_clean_1.fq.gz"
            clean_fq2="${clean_dir}/${sample_name}_clean_2.fq.gz"
            output_dir="${clean_dir}/${sample_name}_output"

            log "处理样本: $sample_name"

            # Step 1: fastp 质量控制
            log "  [fastp] 清洗 reads..."
            fastp -i "$fq1" -I "$fq2" \
                  -o "$clean_fq1" -O "$clean_fq2" \
                  --thread "$threads" --detect_adapter_for_pe \
                  --html "${clean_dir}/${sample_name}_fastp.html" \
                  --json "${clean_dir}/${sample_name}_fastp.json" \
                  || { log "  [错误] fastp 处理失败: $sample_name"; continue; }

            # Step 2: kallisto 定量
            log "  [kallisto] 运行表达量定量..."
            kallisto quant -i "$index_file" -o "$output_dir" \
                "$clean_fq1" "$clean_fq2" -t "$threads" \
                || { log "  [错误] kallisto 运行失败: $sample_name"; continue; }

            log "  [完成] $sample_name"
        else
            log "  [跳过] 未找到 _1/_2.fq.gz: $sample_folder"
        fi
    fi
done

log "所有样本处理完毕"

