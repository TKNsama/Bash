#!/bin/bash

# === Configuration ===
BASE_DIR="${BASE_DIR:-/filer-5/agruppen/PBP/tan}"

# --- 核心配置 ---
DATA_ROOT="$BASE_DIR/wgs/mfs2_Ravi"
WORKSPACE="$BASE_DIR/wgs/mfs2_Ravi/workspace"
REFGENOME="$BASE_DIR/barley/Hv_Morex_v3_hisat2_index/Hv_Morex_v3_pseudomolecules.fa"

SCRIPT_DIR="$WORKSPACE/scripts"
LOG_DIR="$WORKSPACE/logs"
mkdir -p "$SCRIPT_DIR" "$LOG_DIR"
mkdir -p "$WORKSPACE"

# 自动扫描纯数字的样本文件夹
samples=$(find "$DATA_ROOT" -maxdepth 1 -type d -not -path "$DATA_ROOT" | xargs -n 1 basename)

for sample in $samples
do
    # 略过非样本的隐藏文件夹或干扰文件夹
    [[ "$sample" =~ ^[0-9]+$ ]] || continue

    echo "Generating BSA sbatch for: $sample"
    sbatch_file="$SCRIPT_DIR/bsa_${sample}.sbatch"

    cat <<EOT > "$sbatch_file"
#!/bin/bash
#SBATCH --job-name=BSA_${sample}
#SBATCH --auks=yes
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=300G
#SBATCH --output=$LOG_DIR/${sample}_%j.out
#SBATCH --error=$LOG_DIR/${sample}_%j.err
#SBATCH --time=180:11:42

# 1. 加载模块
module load bwa
module load jdk
module load samtools
module load gatk

export APPTAINER_BIND="/filer-5"

# 2. 软件与路径定义
GATK="\${GATK_CONTAINER}"
SAMPLE_DIR="$WORKSPACE/$sample"
mkdir -p "\$SAMPLE_DIR/tmp"

# 严格错误控制：一旦出错立即退出，防止残包糊弄过关
set -e
set -o pipefail

# --- 【动态捕获长文件名】 ---
FQ1=\$(ls $DATA_ROOT/$sample/*_1.fq.gz | head -n 1)
FQ2=\$(ls $DATA_ROOT/$sample/*_2.fq.gz | head -n 1)

# --- 【文件检测环节】 ---
if [ -s "\$SAMPLE_DIR/${sample}_dedup.bam" ]; then
    echo "=========================================================="
    echo "--> [CHECKPOINT] ${sample}_dedup.bam already exists."
    echo "--> Skipping Alignment and Dedup. Proceeding directly to HaplotypeCaller..."
    echo "=========================================================="
else
    echo "Starting BSA Alignment Pipeline for $sample at \$(date)"

    # --- 优化后的 Step A & B & C 整合流 ---
    # 核心逻辑：利用 samtools fixmate -m 给 name-sorted 数据添加错配标签，
    # 接着只进行一次坐标排序(samtools sort)，最后直接 markdup。

    echo "--> Step 1: Running BWA mem & samtools fixmate..."
    # 32线程比对，4线程专门做向后传递
    bwa mem -t 32 -M -R "@RG\tID:$sample\tSM:$sample\tPL:MGI\tLB:lib1" \\
        "$REFGENOME" \\
        "\$FQ1" \\
        "\$FQ2" | \\
        \${SAMTOOLS_CONTAINER} fixmate -@ 4 -m -u - "\$SAMPLE_DIR/${sample}_fixmate.bam"

    echo "--> Step 2: Sorting BAM with strict memory limit..."
    # 32 线程 * 每个线程 5G = 160G 内存，剩余 140G 留给系统缓存，极其安全。
    # -T 指定临时文件写在样本目录下，避免塞爆节点默认的 /tmp
    \${SAMTOOLS_CONTAINER} sort -@ 32 -m 5G \\
        -T "\$SAMPLE_DIR/tmp_sort" \\
        -o "\$SAMPLE_DIR/${sample}_sorted.bam" \\
        "\$SAMPLE_DIR/${sample}_fixmate.bam"

    echo "--> Step 3: Marking duplicates..."
    # 从排好序的文件直接做 markdup，-r 代表移去重复 Reads
    \${SAMTOOLS_CONTAINER} markdup -@ 16 -r \\
        "\$SAMPLE_DIR/${sample}_sorted.bam" \\
        "\$SAMPLE_DIR/${sample}_dedup.bam"

    # 及时清理过程文件，释放上百G的中间磁盘空间
    rm -f "\$SAMPLE_DIR/${sample}_fixmate.bam" "\$SAMPLE_DIR/${sample}_sorted.bam"
fi

# ==========================================================
# 变异呼叫阶段 (BSA DNA-Seq 标准)
# ==========================================================

echo "--> Ensuring CSI index for dedup BAM..."
\${SAMTOOLS_CONTAINER} index -@ 10 -c "\$SAMPLE_DIR/${sample}_dedup.bam"

# --- Step D: HaplotypeCaller ---
echo "--> Running GATK HaplotypeCaller..."
# 给 GATK 分配 48g 内存
\$GATK --java-options "-Xmx48g" HaplotypeCaller \\
    -R "$REFGENOME" \\
    -I "\$SAMPLE_DIR/${sample}_dedup.bam" \\
    -O "\$SAMPLE_DIR/${sample}_raw.vcf" \\
    -L chr1H -L chr2H -L chr3H -L chr4H -L chr5H -L chr6H -L chr7H \\
    --standard-min-confidence-threshold-for-calling 30.0 \\
    --native-pair-hmm-threads 16

# --- Step E: 基础过滤 (GATK DNA 硬过滤标准) ---
echo "--> Running GATK VariantFiltration..."
\$GATK --java-options "-Xmx16g" VariantFiltration \\
    -R "$REFGENOME" \\
    -V "\$SAMPLE_DIR/${sample}_raw.vcf" \\
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \\
    --filter-name "BSA_DNA_HARD_FILTER" \\
    -O "\$SAMPLE_DIR/${sample}_filtered.vcf" \\
    --tmp-dir "\$SAMPLE_DIR/tmp"

# 清理 GATK 的临时小碎片
rm -rf "\$SAMPLE_DIR/tmp"/*

echo "Finished BSA pipeline for $sample at \$(date)"
EOT

    # 提交脚本
    sbatch "$sbatch_file"
done


