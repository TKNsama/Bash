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

# 检查是否成功加载模块
module load blast+
module load bedtools 

# 设置参数
query_file="/filer-5/agruppen/PBP/tan/poaceae_database/query.fa"
output_res="blastp_output.txt"
output_dir="/filer-5/agruppen/PBP/tan/poaceae_database/res"
work_dir="/filer-5/agruppen/PBP/tan/poaceae_database/amino_acid"

cd "$work_dir" || exit 1

# 进行 BLAST
> "${output_dir}/${output_res}"  # 清空或创建最终输出文件
for db_file in "$work_dir"/*.fasta; do
    base_name=$(basename "$db_file")
    echo "Processing $base_name..."

    # 运行 BLAST
    blastp -query "$query_file" -db "$db_file" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq staxid ssciname" \
        -out "${output_dir}/${base_name}_blast_results.txt" -max_target_seqs 15

    # 添加文件名到结果中
    sed -i "1s/^/$base_name\n/" "${output_dir}/${base_name}_blast_results.txt"

    # 拼接结果
    cat "${output_dir}/${base_name}_blast_results.txt" >> "${output_dir}/${output_res}"
done

# 切换到输出目录
cd "$output_dir" || exit 1

# 清理中间结果（如果需要保留中间结果，可以注释掉此行）
rm *_blast_results.txt

# 提取基因名字
awk 'NF && !/^#/{print $2}' "$output_res" | sort -u > name.temp

# 输入参数
output_file="blast_res_seq.fasta"  # 输出文件名
gene_list="name.temp"  # 包含要提取的基因名字的文件

# 检查输出文件是否存在，存在则清空，不存在则创建
> "$output_file"

# 遍历 fasta 文件并提取基因序列
while IFS= read -r gene_name; do
    echo "提取基因 $gene_name 的序列..."
    awk -v gene="^>$gene_name\$" '/^>/{p=0} $0 ~ gene{p=1} p' *.fasta >> "$output_file"
done < "$gene_list"

echo "提取完成。结果保存在 $output_file 中。"

# 清理临时文件
rm name.temp
