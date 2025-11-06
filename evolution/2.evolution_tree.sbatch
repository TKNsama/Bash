#!/bin/bash
 
#SBATCH --job-name=fasttree                                 # Job name, will show up in squeue output
#SBATCH --auks=yes                                      # use Kerberos
#SBATCH --partition=cpu                                 # possible values: cpu, gpu
#SBATCH --ntasks=1                                      # Number of tasks (default=1)
#SBATCH --cpus-per-task=8                               # number of cpus for this task (default=1)
#SBATCH --mem=100G                                       # size of memory (default 10G in cpu and gpu)
#SBATCH --nodes=1                                       # Ensure that all cores are on one machine
#SBATCH --output=job_%j.out                             # File to which standard out will be written
#SBATCH --error=job_%j.err                              # File to which standard err will be written
#SBATCH --mail-type=END                                 # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=tan@ipk-gatersleben.de  # Email to which notifications will be sent 
# load module

module load mafft
module load iqtree
module load seqkit
module load EMBOSS
module load trimAl
# env of pal2nal
PAL2NAL="/filer-5/agruppen/PBP/tan/software/pal2nal.v14/pal2nal.pl"
# the input file from 1.cds_cleaner
fasta="ancient_clade.id.cds.clean.cds.fasta"

# 1) Translate CDS → Protein
transeq -sequence "$fasta" -outseq "${fasta}.pep" -table 1 -frame 1 -clean
sed -i 's/_1$//' ${fasta}.pep
sed -i 's/X$/*/g' ${fasta}.pep
# 2) Align protein (fast, reliable)
mafft --auto "${fasta}.pep" > "${fasta}.pep.aln"

# 3) Back-translate to codon alignment
$PAL2NAL \
    "${fasta}.pep.aln" \
    "$fasta" \
    -output fasta \
    > "${fasta}.codon.fas"

trimal -in "${fasta}.codon.fas" -out "${fasta}.codon.trimmed.fas" -gappyout

# 专注于获得能运行的树，不追求完美bootstrap
iqtree2 -s "${fasta}.codon.trimmed.fas" \
        -st CODON \
        -m MFP \
        -nt 4 \
	-bb 1000 \
        -eps 0.1 \
        -pre ${fasta}.hyphy_tree


