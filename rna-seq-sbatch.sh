raw_data_dir="/filer-5/user/tan/barley_expression_database_rawdata/download.big.ac.cn/gsa2/CRA011642"
processed_data_dir="/filer-5/user/tan/barley_expression_database_rawdata/01fp"
index_dir="/filer-5/user/tan/indexDir"
pair1="f1"
pair2="r2"

echo "#!/bin/bash

#SBATCH --array=1-55
#SBATCH --job-name=rnaSeqStep1Fastp                                 # Job name, will show up in squeue output
#SBATCH --auks=yes                                      # use Kerberos
#SBATCH --partition=cpu                                 # possible values: cpu, gpu
#SBATCH --ntasks=1                                      # Number of tasks (default=1)
#SBATCH --cpus-per-task=4                               # number of cpus for this task (default=1)
#SBATCH --mem=20G                                       # size of memory (default 10G in cpu and gpu)
#SBATCH --nodes=1                                       # Ensure that all cores are on one machine
#SBATCH --output=job_%j.out                             # File to which standard out will be written
#SBATCH --error=job_%j.err                              # File to which standard err will be written
#SBATCH --mail-type=END                                 # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=tan@ipk-gatersleben.de  # Email to which notifications will be sent

# don't know why, but it should be inputed
source /etc/profile
# Load necessary modules
module load fastp
module load kallisto/0.46.1
" > rnaSeq.sbatch

# get name
for sample_folder in "$raw_data_dir"/*; do
    # Check if it's a directory
    if [ -d "$sample_folder" ]; then
        # Get the sample name based on the folder name
        sample_name=$(basename "$sample_folder")
        echo "fastp -i \"$sample_folder/${sample_name}_${pair1}.fq.gz\" -I \"$sample_folder/${sample_name}_${pair2}.fq.gz\" -o \"$processed_data_dir/kallisto/${sample_name}_clean_${pair1}.fq.gz\" -O \"$processed_data_dir/kallisto/${sample_name}_clean_${pair2}.fq.gz\" --thread 8" >> rnaSeq.sbatch
    fi
done

for sample_folder in "$raw_data_dir"/*; do
    # Check if it's a directory
    if [ -d "$sample_folder" ]; then
        # Get the sample name based on the folder name
        sample_name=$(basename "$sample_folder")
        echo "kallisto quant -i \"$index_dir/190823_Barley_Morex_V2_gene_annotation_PGSB.all.cds.annotations.idx\" -o \"$processed_data_dir/kallisto/${sample_name}_output\" \"$processed_data_dir/kallisto/${sample_name}_clean_${pair1}.fq.gz\" \"$processed_data_dir/kallisto/${sample_name}_clean_${pair2}.fq.gz\" --paired -t 8" >> rnaSeq.sbatch
    fi
done
