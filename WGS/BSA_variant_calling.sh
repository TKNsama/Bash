#!/bin/bash

#SBATCH --job-name=index                                 # Job name, will show up in squeue output
#SBATCH --auks=yes                                      # use Kerberos
#SBATCH --partition=cpu                                 # possible values: cpu, gpu
#SBATCH --ntasks=1                                      # Number of tasks (default=1)
#SBATCH --cpus-per-task=50                               # number of cpus for this task (default=1)
#SBATCH --mem=300G                                       # size of memory (default 10G in cpu and gpu)
#SBATCH --nodes=1                                       # Ensure that all cores are on one machine
#SBATCH --output=job_%j.out                             # File to which standard out will be written
#SBATCH --error=job_%j.err                              # File to which standard err will be written
#SBATCH --mail-type=END                                 # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=tan@ipk-gatersleben.de  # Email to which notifications will be sent
#SBATCH --time=168:00:00

# Load necessary modules

module load bwa
module load samtools
module load sambamba
module load vcftools
module load gatk

# Define workspace directory
workspace="/filer-5/agruppen/PBP/tan/wgs/workspace/bowman"

# Define reference genome and sample list
refgenome="/filer-5/agruppen/PBP/tan/indexDir/barley_index/220830_Bowman_pseudomolecules_and_unplaced_contigs_CPclean.fasta"
samples=("sample1" "sample2" "sample3" "sample4") # Replace with your actual sample names
name="KENANJAN2025"
log_file="$workspace/pipeline.log"
exec > >(tee -i "$log_file") 2>&1

# Step 1: Prepare reference genome
if [[ ! -f "${refgenome%.fasta}.dict" ]]; then
  java -jar /opt/Bio/picard-tools/2.21.9/picard.jar CreateSequenceDictionary \
      R="$refgenome" \
      O="${refgenome%.fasta}.dict"
fi
if [[ ! -f "${refgenome}.fai" ]]; then
  samtools faidx "$refgenome"
fi

# Step 2 to Step 14: Process each sample in the list
for sample in "${samples[@]}"
do
  echo "Processing $sample"

  # Create output directories
  sample_dir="$workspace/$sample"
  mkdir -p "$sample_dir"
  bwa_dir="$sample_dir/bwa"
  picard_dir="$sample_dir/picard"
  gatk_dir="$sample_dir/gatk"
  mkdir -p "$bwa_dir" "$picard_dir" "$gatk_dir"

  # Step 2: Alignment of reads to the reference genome (if BAM file doesn't exist)
  bam_file="$bwa_dir/${sample}_sorted.bam"
  if [[ -f "$bam_file" ]]; then
    echo "BAM file for $sample already exists. Skipping alignment and BAM conversion."
  else
    bwa mem -t 48 -M -R "@RG\tPU:"${name}"\tID:"${sample}"\tSM:"${sample}"\tLB:WXS\tPL:ILLUMINA" "$refgenome" "/filer-5/agruppen/PBP/tan/wgs/workspace/${sample}_1.fq.gz" "/filer-5/agruppen/PBP/tan/wgs/workspace/${sample}_2.fq.gz" | \
    sambamba view -S -t 48 -f bam /dev/stdin | \
    sambamba sort -t 48 -n --tmpdir="$workspace/tmp" -o "$bwa_dir/${sample}_sorted.bam" /dev/stdin
  fi

  gatk MarkDuplicates -I "$bwa_dir/${sample}_sorted.bam" -O "$picard_dir/${sample}_marked.bam" -M "$picard_dir/${sample}_duplication_metrics.txt"

  # Step 4: Check the output of the alignment using SAMtools flagstat
  samtools flagstat "$bwa_dir/${sample}_sorted.bam" > "$bwa_dir/${sample}_alignment_metrics.txt"

  # Step 6: Mark duplicates in the BAM file using Sambamba
  sambamba sort -t 48 "$picard_dir/${sample}_marked.bam" -o "$picard_dir/${sample}_sorted_marked.bam"
  sambamba index -c "$picard_dir/${sample}_sorted_marked.bam"

  # Step 7: Call variants using GATK HaplotypeCaller
  gatk HaplotypeCaller \
      -R "$refgenome" \
      -I "$picard_dir/${sample}_sorted_marked.bam" \
      -O "$gatk_dir/raw_${sample}.vcf"

  # Step 8: Extract SNPs from the raw VCF
  gatk SelectVariants \
      -R "$refgenome" \
      -V "$gatk_dir/raw_${sample}.vcf" \
      --select-type-to-include SNP \
      -O "$gatk_dir/${sample}_raw_snps.vcf"

  # Step 9: Extract Indels from the raw VCF
  gatk SelectVariants \
      -R "$refgenome" \
      -V "$gatk_dir/raw_${sample}.vcf" \
      --select-type-to-include INDEL \
      -O "$gatk_dir/${sample}_raw_indels.vcf"

  # Step 10: Filter the SNPs based on quality metrics
  gatk VariantFiltration \
      -R "$refgenome" \
      -V "$gatk_dir/${sample}_raw_snps.vcf" \
      --filter-expression "QD < 2.0 || FS > 60.0 || SOR > 10.0 || MQ < 40.0" \
      --filter-name "MY_filter" \
      -O "$gatk_dir/filtered_${sample}_snps.vcf"

  echo "$sample processing complete."
done

# Step 11: Merge all VCF files using vcftools
merged_dir="$workspace/merged"
mkdir -p "$merged_dir"
vcf-merge "$workspace"/*/*/*.vcf > "$merged_dir/final_merged.vcf"

echo "VCF files merged into $merged_dir/final_merged.vcf"
