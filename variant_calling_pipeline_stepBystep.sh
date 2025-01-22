###First try ipyrad.
###then try following method
###we will compare the results
###please replace the s1 with correct name of one of the samples
###please replace the reference.fasta with the correct fererence file name (use the masked reference)


# Step 1: Generate the index of the reference genome
module load bwa
bwa index -a bwtsw reference.fasta

# Step 2: Alignment of reads to the reference genome

bwa mem -M reference.fasta s1.fastq > s1.sam

# Step 3: Check the output of the alignment
module load samtools
samtools flagstat s1.sam > alignment_metrics.txt

# Step 4: Convert SAM to BAM and sort by genomic position
module load jdk
java -jar /opt/Bio/picard-tools/2.21.9/picard.jar SortSam \
    INPUT=s1.sam \
    OUTPUT=sorted_s1.bam \
    SORT_ORDER=coordinate

# Step 5: Add Read Groups (RG) to BAM file
java -jar /opt/Bio/picard-tools/2.21.9/picard.jar AddOrReplaceReadGroups \
    INPUT=sorted_s1.bam \
    OUTPUT=sorted_s1_RG.bam \
    RGPL=illumina \
    RGPU=D109LACXX \
    RGLB=Lib1 \
    RGID=s1 \
    RGSM=s1 \
    VALIDATION_STRINGENCY=LENIENT

# Step 6: Mark duplicates in the BAM file
java -jar /opt/Bio/picard-tools/2.21.9/picard.jar MarkDuplicates \
    INPUT=sorted_s1_RG.bam \
    OUTPUT=dedup_sorted_s1.bam \
    METRICS_FILE=metrics.txt

# Step 7: Create a sequence dictionary for the reference genome
java -jar /opt/Bio/picard-tools/2.21.9/picard.jar CreateSequenceDictionary \
    R=reference.fasta \
    O=reference.dict

# Step 8: Generate index files for the reference genome and the BAM file
samtools faidx reference.fasta
samtools index dedup_sorted_s1.bam

# Step 9: Call variants using GATK HaplotypeCaller
java -Xmx6g -jar /opt/Bio/gatk/4.1.9.0/gatk-package-4.1.9.0-local.jar HaplotypeCaller \
    -R reference.fasta \
    -I dedup_sorted_s1.bam \
    -O raw_s1.vcf

# Step 10: Generate an index from a known-sites VCF file (if using BQSR)
# This is required for base quality score recalibration.
java -Xmx6g -jar /opt/Bio/gatk/4.1.9.0/gatk-package-4.1.9.0-local.jar IndexFeatureFile \
    -I wheat.vcf

# Step 11: Base Quality Score Recalibration (BQSR)
java -Xmx6g -jar /opt/Bio/gatk/4.1.9.0/gatk-package-4.1.9.0-local.jar BaseRecalibrator \
    -I dedup_sorted_s1.bam \
    -R reference.fasta \
    -O recal_data.table \
    --known-sites wheat.vcf

# Step 12: Extract SNPs from the raw VCF
java -Xmx6g -jar /opt/Bio/gatk/4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants \
    -R reference.fasta \
    -V raw_s1.vcf \
    --select-type-to-include SNP \
    -O s1_raw_snps.vcf

# Step 13: Extract Indels from the raw VCF
java -Xmx6g -jar /opt/Bio/gatk/4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants \
    -R reference.fasta \
    -V raw_s1.vcf \
    --select-type-to-include INDEL \
    -O s1_raw_indels.vcf

# Step 14: Filter the SNPs based on quality metrics
java -Xmx6g -jar /opt/Bio/gatk/4.1.9.0/gatk-package-4.1.9.0-local.jar VariantFiltration \
    -R reference.fasta \
    -V s1_raw_snps.vcf \
    --filter-expression "QD < 2.0 || FS > 60.0 || SOR > 10.0 || MQ < 40.0" \
    --filter-name "MY_filter" \
    -O filtered_snps.vcf
#repat these steps for all sample:
you can do it using a for loop. But first run you code one by one for one sample. when there is no error, run the rest using a for loop

#then merge all vcf files using vcftools:
vcf-merge *.vcf > final_merged.vcf
