#!/bin/bash

#SBATCH --job-name=gwas1                          # Job name
#SBATCH --auks=yes                                     # Use Kerberos
#SBATCH --partition=cpu                                # Partition type
#SBATCH --ntasks=1                                     # Number of tasks
#SBATCH --cpus-per-task=16                            # Number of CPUs per task
#SBATCH --mem=100G                                      # Memory allocation
#SBATCH --nodes=1                                      # Ensure all cores are on one machine
#SBATCH --output=job_%j.out                         # Standard output log
#SBATCH --error=job_%j.err                          # Standard error log
#SBATCH --time=200:23:37

module load plink
module load gemma

# ================== Configuration ==================
VCF="gwas.vcf.gz"          # 输入 VCF
PREFIX="gwas"              # PLINK / GEMMA 前缀（不要带扩展名）
N_PC=10                    # PCA 数量

PLINK=plink
GEMMA=gemma

# ================== Step 1: VCF → PLINK ==================
echo "Step 1: Converting VCF to PLINK format..."
$PLINK \
  --vcf "$VCF" \
  --geno 0.1 \
  --maf 0.1 \
  --make-bed \
  --double-id \
  --allow-extra-chr \
  --out "$PREFIX"

