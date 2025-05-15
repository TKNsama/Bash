#!/bin/bash

module load seqkit

# name 文件是基因ID，从blast结果中挑选/提取出来即可

for i in Barley/*.fasta; do 
	cat $i | seqkit grep -f name >> res.fasta
	echo "$i is done"
 done




