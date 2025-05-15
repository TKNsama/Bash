#!/bin/bash
 
#SBATCH --job-name=rnaSeqStep1Fastp                                 # Job name, will show up in squeue output
#SBATCH --auks=yes                                      # use Kerberos
#SBATCH --partition=cpu                                 # possible values: cpu, gpu
#SBATCH --ntasks=1                                      # Number of tasks (default=1)
#SBATCH --cpus-per-task=4                               # number of cpus for this task (default=1)
#SBATCH --mem=20G                                       # size of memory (default 10G in cpu and gpu)
#SBATCH --nodes=1                                       # Ensure that all cores are on one machine
#SBATCH --output=job_%j.out                             # File to which standard out will be written
#SBATCH --error=job_%j.err                              # File to which standard err will be written
#SBATCH --mail-type=END                                 # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=tan@ipk-gatersleben.de  # Email to which notifications will be sent 
module load blast+
module load bedtools


#!/bin/bash

# ===== 用户自定义配置 =====
query_file="/filer-5/agruppen/PBP/tan/database/relative_species_pep/query.fa"
output_res="blast_output.txt"
output_dir="/filer-5/agruppen/PBP/tan/database/relative_species_pep/res"
work_dir="/filer-5/agruppen/PBP/tan/database/relative_species_pep" # blastdb位置
max_hits=5
threads=4  # 可修改线程数

# ===== 环境准备与验证 =====
mkdir -p "$output_dir"
cd "$work_dir" || { echo "无法进入工作目录 $work_dir"; exit 1; }

# 清空旧输出文件
: > "${output_dir}/${output_res}"  # 清空内容但不删除文件

# ===== 批量 BLAST 比对 =====
echo ">>> 开始 BLASTP 比对任务: $(date)"
for db_file in *.fasta *.fa *.faa; do
    [[ -e "$db_file" ]] || continue  # 若无匹配文件则跳过

    base_name=$(basename "$db_file")
    db_name="${base_name%.*}"  # 去掉扩展名
    result_file="${output_dir}/${db_name}_blast_results.tmp.txt"

    echo ">> 比对 $query_file vs $db_file"

    blastp -query "$query_file" \
           -db "$db_file" \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq staxid ssciname" \
           -max_target_seqs "$max_hits" \
           -num_threads "$threads" \
           -out "$result_file"

    if [[ -s "$result_file" ]]; then
        # 插入数据库文件名作为前缀列
        sed -i "1s/^/${db_name}\t/" "$result_file"

        # 取前5行追加到主文件
        head -n 5 "$result_file" >> "${output_dir}/${output_res}"
    else
        echo "结果为空: $db_file" >&2
    fi
done

# ===== 清理中间文件 =====
rm -f "${output_dir}"/*_blast_results.tmp.txt

echo ">>> BLASTP 比对完成: $(date)"
echo "输出文件: ${output_dir}/${output_res}"


