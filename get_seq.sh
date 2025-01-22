#!/bin/bash

module load seqkit

for i in amino_acid/*.fasta; do 
	cat $i | seqkit grep -f name >> res.fasta
	echo "$i is done"
 done
