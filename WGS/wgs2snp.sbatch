#!/bin/bash
#SBATCH --job-name=wgs_mapping_snp
#SBATCH --auks=yes
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=tan@ipk-gatersleben.de

module load bwa
module load samtools
module load jdk
module load vcftools

# work space (sample location)
workspace="/filer-5/agruppen/PBP/tan/wgs/mfs1_F2_mapping"

# ref genome
refgenome="/filer-5/agruppen/PBP/tan/indexDir/barley_index/220809_Foma_pseudomolecules_and_unplaced_contigs_CPclean.fasta"

log_file="$workspace/pipeline.log"
exec > >(tee -i "$log_file") 2>&1

# dict creation of ref genome
if [[ ! -f "${refgenome%.fasta}.dict" ]]; then
  java -jar /opt/Bio/picard-tools/2.21.9/picard.jar CreateSequenceDictionary \
      R="$refgenome" \
      O="${refgenome%.fasta}.dict"
fi
if [[ ! -f "${refgenome}.fai" ]]; then
  samtools faidx "$refgenome"
fi

# processing
for sample_dir in "$workspace"/*; do
  [[ -d "$sample_dir" ]] || continue
  sample_id=$(basename "$sample_dir")

  # check fastaq
  fq_files=("$sample_dir"/2/*.fastq.gz)
  if [[ ${#fq_files[@]} -eq 0 ]]; then
    echo "Warning: No FASTQ found for $sample_id"
    continue
  fi

  # R1 or R1&R2
  fq1="${fq_files[0]}"
  fq2=""
  for f in "${fq_files[@]}"; do
    [[ "$f" == *"_R2_"* ]] && fq2="$f"
  done

  echo "Processing sample: $sample_id"
  bwa_out="$sample_dir/bwa"
  picard_out="$sample_dir/picard"
  gatk_out="$sample_dir/gatk"
  mkdir -p "$bwa_out" "$picard_out" "$gatk_out"

  # alignment
  if [[ -n "$fq2" ]]; then
    bwa mem -t 48 -M "$refgenome" "$fq1" "$fq2" > "$bwa_out/$sample_id.sam"
  else
    bwa mem -t 48 -M "$refgenome" "$fq1" > "$bwa_out/$sample_id.sam"
  fi

  samtools flagstat "$bwa_out/$sample_id.sam" > "$bwa_out/${sample_id}_alignment_metrics.txt"

  java -jar /opt/Bio/picard-tools/2.21.9/picard.jar SortSam \
      INPUT="$bwa_out/$sample_id.sam" \
      OUTPUT="$picard_out/sorted_${sample_id}.bam" \
      SORT_ORDER=coordinate

  java -jar /opt/Bio/picard-tools/2.21.9/picard.jar AddOrReplaceReadGroups \
      INPUT="$picard_out/sorted_${sample_id}.bam" \
      OUTPUT="$picard_out/sorted_${sample_id}_RG.bam" \
      RGPL=MGI \
      RGPU=D109LACXX \
      RGLB=Lib1 \
      RGID="$sample_id" \
      RGSM="$sample_id" \
      VALIDATION_STRINGENCY=LENIENT

  java -jar /opt/Bio/picard-tools/2.21.9/picard.jar MarkDuplicates \
      INPUT="$picard_out/sorted_${sample_id}_RG.bam" \
      OUTPUT="$picard_out/dedup_sorted_${sample_id}.bam" \
      METRICS_FILE="$picard_out/${sample_id}_metrics.txt"
  
  samtools index -c "$picard_out/dedup_sorted_${sample_id}.bam"

  java -Xmx16g -jar /opt/Bio/gatk/4.1.9.0/bin/gatk-package-4.1.9.0-local.jar HaplotypeCaller \
      -R "$refgenome" \
      -I "$picard_out/dedup_sorted_${sample_id}.bam" \
      -O "$gatk_out/raw_${sample_id}.vcf"

  java -Xmx16g -jar /opt/Bio/gatk/4.1.9.0/bin/gatk-package-4.1.9.0-local.jar SelectVariants \
      -R "$refgenome" \
      -V "$gatk_out/raw_${sample_id}.vcf" \
      --select-type-to-include SNP \
      -O "$gatk_out/${sample_id}_raw_snps.vcf"

  java -Xmx16g -jar /opt/Bio/gatk/4.1.9.0/bin/gatk-package-4.1.9.0-local.jar SelectVariants \
      -R "$refgenome" \
      -V "$gatk_out/raw_${sample_id}.vcf" \
      --select-type-to-include INDEL \
      -O "$gatk_out/${sample_id}_raw_indels.vcf"

  java -Xmx16g -jar /opt/Bio/gatk/4.1.9.0/bin/gatk-package-4.1.9.0-local.jar VariantFiltration \
      -R "$refgenome" \
      -V "$gatk_out/${sample_id}_raw_snps.vcf" \
      --filter-expression "QD < 2.0 || FS > 60.0 || SOR > 10.0 || MQ < 40.0" \
      --filter-name "MY_filter" \
      -O "$gatk_out/filtered_${sample_id}_snps.vcf"

  echo "$sample_id processing complete."
done

