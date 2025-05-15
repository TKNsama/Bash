#!/bin/bash


# Loop to generate and submit jobs for each index
for i in {1..21}; do
    script_name="gemma_job_$i.sh"
    echo "#!/bin/bash" > $script_name
    echo "#SBATCH --job-name=gemma_$i" >> $script_name
    echo "#SBATCH --auks=yes" >> $script_name
    echo "#SBATCH --output=job_%j.out" >> $script_name
    echo "#SBATCH --error=job_%j.err" >> $script_name
    echo "#SBATCH --partition=cpu" >> $script_name
    echo "#SBATCH --cpus-per-task=24" >> $script_name
    echo "#SBATCH --mem=100G" >> $script_name
    echo "#SBATCH --nodes=1" >> $script_name
    echo "" >> $script_name
    echo "date" >> $script_name
    echo "" >> $script_name
    echo "module load gemma" >> $script_name
    echo "" >> $script_name
    echo "gemma -bfile grains_248_filtered.vcf \\" >> $script_name
    echo "      -k output/kinship.cXX.txt \\" >> $script_name
    echo "      -n $i \\" >> $script_name
    echo "      -lmm 1 \\" >> $script_name
    echo "      -c snp_covariates.covar \\" >> $script_name   # optional
    echo "      -o ${i}_res" >> $script_name
    echo "" >> $script_name
    echo "date" >> $script_name

    sbatch $script_name
    echo "$i is submitted"
done
