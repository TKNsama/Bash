#!/bin/bash
#SBATCH --job-name=gwas1 # Job name
#SBATCH --auks=yes # Use Kerberos
#SBATCH --partition=cpu # Partition type
#SBATCH --ntasks=1 # Number of tasks
#SBATCH --cpus-per-task=16 # Number of CPUs per task
#SBATCH --mem=100G # Memory allocation
#SBATCH --nodes=1 # Ensure all cores are on one machine
#SBATCH --output=job_%j.out # Standard output log
#SBATCH --error=job_%j.err # Standard error log
#SBATCH --time=200:23:37

module load plink
module load gemma

# ================== Configuration ==================
VCF="gwas.vcf.gz" # 输入 VCF
PREFIX="gwas.vcf.gz" # PLINK / GEMMA 前缀（不要带扩展名）
N_PC=10 # PCA 数量
PLINK=plink
GEMMA=gemma
# ================== Step 2: PCA ==================
echo "Step 2: Running PCA (${N_PC} PCs)..."
$PLINK 
  --bfile "$PREFIX" 
  --pca $N_PC 
  --allow-extra-chr 
  --out "${PREFIX}_PCA"
# ================== Step 3: Format PCA covariates ==================
echo "Step 3: Formatting PCA covariates for GEMMA..."
# PLINK eigenvec: FID IID PC1 PC2 ...
# GEMMA covar: PC1 PC2 ...
awk -v npc=$N_PC '{
  for (i=3; i<3+npc; i++) {
    printf "%s%s", $i, (i==2+npc ? ORS : OFS)
  }
}' "${PREFIX}_PCA.eigenvec" > pca_covariates.covar
# ================== Step 4: Kinship matrix ==================
echo "Step 4: Calculating kinship matrix (centered)..."
$GEMMA 
  -bfile "$PREFIX" 
  -gk 1 
  -o kinship
# ================== Step 5: Final checks ==================
echo "Step 5: Checking output files..."
for f in 
  "${PREFIX}.bed" 
  "${PREFIX}.bim" 
  "${PREFIX}.fam" 
  "output/kinship.cXX.txt" 
  "pca_covariates.covar"
do
  [[ -f "$f" ]] || { echo "ERROR: Missing $f"; exit 1; }
done
echo "=============================================="
echo " All GEMMA input files successfully generated "
echo "=============================================="
echo "PLINK prefix : $PREFIX"
echo "Kinship : output/kinship.cXX.txt"
echo "Covariates : pca_covariates.covar"
