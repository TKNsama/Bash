#!/bin/bash
 
#SBATCH --job-name=rnaSeqStep1Fastp                                 # Job name, will show up in squeue output
#SBATCH --auks=yes                                      # use Kerberos
#SBATCH --partition=cpu                                 # possible values: cpu, gpu
#SBATCH --ntasks=1                                      # Number of tasks (default=1)
#SBATCH --cpus-per-task=32                               # number of cpus for this task (default=1)
#SBATCH --mem=50G                                       # size of memory (default 10G in cpu and gpu)
#SBATCH --nodes=1                                       # Ensure that all cores are on one machine
#SBATCH --output=job_%j.out                             # File to which standard out will be written
#SBATCH --error=job_%j.err                              # File to which standard err will be written
#SBATCH --mail-type=END                                 # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=tan@ipk-gatersleben.de  # Email to which notifications will be sent 

# load module

module load FastTree
module load mafft
module load iqtree

fasta="ARF_pan_genome.fasta" # 输入fasta文件名

mafft --auto $fasta > ${fasta}.fas
#FastTree Eukaryota_tar_aa.fas > res_tree.nwk
iqtree2 -s ${fasta}.fas -nt AUTO -T 32

