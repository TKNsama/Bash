#!/bin/bash

#SBATCH --job-name=${prefix}                          # Job name
#SBATCH --auks=yes                                     # Use Kerberos
#SBATCH --partition=cpu                                # Partition type
#SBATCH --ntasks=1                                     # Number of tasks
#SBATCH --cpus-per-task=32                            # Number of CPUs per task
#SBATCH --mem=100G                                      # Memory allocation
#SBATCH --nodes=1                                      # Ensure all cores are on one machine
#SBATCH --output=job_%j.out                         # Standard output log
#SBATCH --error=job_%j.err                          # Standard error log
#SBATCH --exclude=slurm-03,slurm-20                    # Exclude nodes

# ======= Input Files =======
region_file="regions.txt"           # 每行格式: chr start end
vcf_file="filtered_grains_248.vcf.gz"  # 使用过滤后的VCF
tmp_vcf="tmp_covariates.vcf.gz"

# ======= Step 1: Extract SNPs from specified regions =======
echo "Extracting SNPs from specified regions..."
bcftools view -R "$region_file" "$vcf_file" -Oz -o "$tmp_vcf"
bcftools index "$tmp_vcf"

# ======= Step 2: Convert to PLINK format =======
plink_prefix="snp_covariates"
plink --vcf "$tmp_vcf" --make-bed --double-id --out "$plink_prefix"

# ======= Step 3: Convert to .raw format (additive coding) =======
plink --bfile "$plink_prefix" --recode A --out "$plink_prefix"

# ======= Step 4: Format for GEMMA covariate file =======
# Remove header and non-genotype columns, keep only numeric genotype matrix
awk 'NR==1 {for (i=1; i<=NF; i++) colname[i]=$i; next} 
     {printf $1; for (i=7; i<=NF; i++) printf " %s", $i; print ""}' "${plink_prefix}.raw" > snp_covariates.covar

echo "SNP covariate file generated: snp_covariates.covar"
