#!/bin/bash

#SBATCH --job-name=rnaSeqStep1Fastp                                 # Job name, will show up in squeue output
#SBATCH --auks=yes                                      # use Kerberos
#SBATCH --partition=cpu                                 # possible values: cpu, gpu
#SBATCH --ntasks=1                                      # Number of tasks (default=1)
#SBATCH --cpus-per-task=4                               # number of cpus for this task (default=1)
#SBATCH --mem=50G                                       # size of memory (default 10G in cpu and gpu)
#SBATCH --nodes=1                                       # Ensure that all cores are on one machine
#SBATCH --output=job_%j.out                             # File to which standard out will be written
#SBATCH --error=job_%j.err                              # File to which standard err will be written
#SBATCH --mail-type=END                                 # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=tan@ipk-gatersleben.de  # Email to which notifications will be sent

#!/bin/bash

# 加载 BLAST+ 模块
module load blast+

# 参数定义
dbtype="nucl"
input_dir="/filer-5/user/tan/pangenome/selected_library_files/Barley/Barley Pangenome v2 2023/Assemblies/compressed assemblies (gzipped)"

# 检查输入目录是否存在
if [ ! -d "$input_dir" ]; then
    echo "错误：目录不存在：$input_dir"
    exit 1
fi

cd "$input_dir" || { echo "无法进入目录：$input_dir"; exit 1; }

# 日志文件
log_file="makeblastdb_log.txt"
echo "=== BLAST 数据库构建日志 ===" > "$log_file"
echo "开始时间：$(date)" >> "$log_file"
echo >> "$log_file"

# 构建数据库
for file in *.fasta *.fa *.fna; do
    [ -e "$file" ] || continue  # 如果没有匹配文件则跳过

    base_name="${file%.*}"
    echo "正在处理：$file ..." | tee -a "$log_file"

    makeblastdb -in "$file" -dbtype "$dbtype" -parse_seqids >> "$log_file" 2>&1

    if [ $? -eq 0 ]; then
        echo "成功构建数据库：$base_name" >> "$log_file"
    else
        echo "失败：$base_name" >> "$log_file"
    fi

    echo >> "$log_file"
done

echo "完成时间：$(date)" >> "$log_file"
echo "日志保存在：$log_file"




