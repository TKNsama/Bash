#!/bin/bash

# 设置处理器数量（proc 文件路径）
PROC=proc

for cds_file in *_cds.fasta; do
    # 获取前缀
    prefix=${cds_file%_cds.fasta}
    
    # 对应的蛋白文件和 homologs 文件
    pep_file="${prefix}_pep.fasta"
    homolog_file="${prefix}.txt"
    
    # 检查这两个文件是否都存在
    if [[ ! -f "$pep_file" ]]; then
        echo "Missing protein file: $pep_file"
        continue
    fi
    
    if [[ ! -f "$homolog_file" ]]; then
        echo "Missing homologs file: $homolog_file"
        continue
    fi

    # 输出目录
    output_dir="output_${prefix}"
    
    echo "Running ParaAT for $prefix..."

    ./ParaAT.pl -h "$homolog_file" \
              -n "$cds_file" \
              -a "$pep_file" \
              -p "$PROC" \
              -o "$output_dir" \
              -m mafft \
              -f axt 
done

