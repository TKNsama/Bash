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

# load module

module load hmmer

hmmfile="PF06507.hmm"
output_dir="ARF"  # 输出结果目录

# 创建输出结果目录
mkdir -p $output_dir

# 遍历所有 .fa 文件
for i in Barley/*.fasta; do
    # 提取文件名（不包括路径和扩展名）
    filename=$(basename -- "$i")
    filename_no_ext="${filename%.*}"

    # 执行 hmmsearch
    hmmsearch --cut_tc --domtblout "${output_dir}/${filename_no_ext}.out" "$hmmfile" "$i"

    # 如果想要同时输出详细信息到终端，可以取消下一行的注释
    # echo "Finished hmmsearch for ${filename}"

done

# 将所有输出文件合并为 hmm.res
cat "${output_dir}"/*.out > hmm.res

# 根据阈值过滤结果，生成 id.txt
grep -v "^#" hmm.res | awk '($7 + 0) < 1E-20' | cut -f1 -d " " | sort -u > id.txt


# 输出完成消息
echo "Analysis completed. Results saved in id.txt."

