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

# ======= Load modules =======
module load bcftools
module load plink
module load gemma

# ======= Configuration =======
vcf_file="gwas.vcf.gz"        # 输入的VCF文件
sample_file="gwas.id"     # 包含样本名的文本文件（每行一个ID）

# ======= Step 1: Filter VCF by samples =======
echo "Filtering VCF file for selected samples..."
filtered_vcf="filtered_${vcf_file}"
bcftools view -S "$sample_file" "$vcf_file" -Oz -o "$filtered_vcf"
bcftools index "$filtered_vcf"

# ======= Step 2: Convert VCF to PLINK format =======
echo "Converting VCF to PLINK format..."
plink_prefix="${filtered_vcf%.vcf.gz}"  # 去掉.vcf.gz后缀
plink --vcf "$filtered_vcf" --make-bed --double-id --allow-extra-chr --out "$plink_prefix"

# ======= Step 3: PCA Analysis =======
echo "Running PCA..."
plink --bfile "$plink_prefix" --pca 10 --out "${plink_prefix}_PC"

# ======= Step 4: Kinship Matrix =======
echo "Calculating kinship matrix using GEMMA..."
gemma -bfile "$plink_prefix" -gk 1 -o kinship_filtered

echo "All steps completed."
