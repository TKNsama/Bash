#!/bin/bash

#SBATCH --job-name=gwas2                          # Job name
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
# ===== Configuration =====

VCF="filtered_gwas.vcf.gz"   # 输入 VCF
OUT="gwas"          # 与 gemma -bfile 保持一致
N_PC=10                             # PCA covariate 数量

PLINK=plink
GEMMA=gemma
# ===== Step 0: Check =====
echo "Checking input files..."
[[ -f "$VCF" ]] || { echo "ERROR: VCF not found"; exit 1; }
[[ -f "$VCF.tbi" ]] || { echo "ERROR: VCF index not found"; exit 1; }
# ===== Step 1: VCF → PLINK =====
echo "Converting VCF to PLINK (wheat chromosomes)..."
$PLINK \
  --vcf "$VCF" \
  --make-bed \
  --double-id \
  --allow-extra-chr \
  --out "$OUT"
# ===== Step 2: PCA → covariates =====
echo "Running PCA..."
$PLINK \
  --bfile "$OUT" \
  --pca $N_PC \
  --out "${OUT}_PCA"

echo "Formatting PCA covariates for GEMMA..."
# eigenvec: FID IID PC1 PC2 ...
awk -v npc=$N_PC '{
  for (i=3; i<3+npc; i++) {
    printf "%s%s", $i, (i==2+npc ? ORS : OFS)
  }
}' "${OUT}_PCA.eigenvec" > snp_covariates.covar
# ===== Step 3: Kinship =====
echo "Calculating kinship matrix..."
$GEMMA \
  -bfile "$OUT" \
  -gk 1 \
  -o kinship
# ===== Step 4: Final check =====
echo "Checking outputs..."

for f in \
  "${OUT}.bed" \
  "${OUT}.bim" \
  "${OUT}.fam" \
  "output/kinship.cXX.txt" \
  "snp_covariates.covar"
do
  [[ -f "$f" ]] || { echo "ERROR: Missing $f"; exit 1; }
done

echo "================================="
echo "All GEMMA input files generated."
echo "You can now submit GWAS jobs."
echo "================================="

