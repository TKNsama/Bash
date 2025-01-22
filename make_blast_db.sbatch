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

module load blast+

dbtype="prot" # prot or nucl
input_dir="/filer-5/agruppen/PBP/tan/database/grasses_genome/species_uniq"

cd "$input_dir"

# 注意：在 Bash 中使用数组时应使用括号 ()
for file in *.fasta; do
    base_name="${file%.fasta}"
    makeblastdb -in "$file" -dbtype "$dbtype" -parse_seqids
done



