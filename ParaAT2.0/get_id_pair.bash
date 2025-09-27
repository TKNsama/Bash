#!/bin/bash

# 遍历所有 *_cds.fasta 文件
for fasta in *_cds.fasta; do
    # 提取前缀作为输出文件名
    base="${fasta%_cds.fasta}"
    outfile="${base}.txt"

    # 提取所有序列名（去掉前缀的 ">"）
    seqs=($(grep "^>" "$fasta" | sed 's/^>//'))

    # 如果序列少于2条，就跳过
    if [[ ${#seqs[@]} -lt 2 ]]; then
        echo "Skipping $fasta (less than 2 sequences)"
        continue
    fi

    # 第一个序列作为 anchor
    anchor="${seqs[0]}"

    # 生成配对文件
    : > "$outfile"  # 清空/新建文件
    for ((i=1; i<${#seqs[@]}; i++)); do
        echo -e "${anchor}\t${seqs[$i]}" >> "$outfile"
    done

    echo "Generated: $outfile"
done

