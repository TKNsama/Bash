#!/bin/bash

# 遍历所有以 output_N0 开头的文件夹
for group in output_*
do
    # 创建对应的输出目录
    mkdir -p kaks/$group

    # 遍历该组下所有 .axt 文件
    for file in $group/*.axt
    do
        # 取得文件名
        fname=$(basename $file)
        # 执行 KaKs_Calculator
        KaKs_Calculator -i $file -o kaks/$group/${fname}.kaks -m YN
    done
done
