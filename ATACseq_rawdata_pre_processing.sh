#!/bin/bash
# ============================================================================
# ChIP-seq pipeline: QC -> trimming -> alignment -> dedup -> peak calling
# ============================================================================
set -euo pipefail

# ------------------------- path -------------------------
RAW_DIR="/storage2/vitoLab/xiongzy12295001/ATACseq/test/K10"  #"/RawData"                          # 原始 fastq 所在目录，按 ${sample}/${sample}_1.fq.gz 组织
WORK_DIR="/storage2/vitoLab/xiongzy12295001/ATACseq/test/work" #"/workspace"                       # 输出目录

FASTQC_BIN="/soft/bio/fastqc-0.11.9/fastqc" #"/path/fastqc"
FASTQ_SCREEN_BIN="/storage/shihongjunLab/liulifeng/workflow/soft/FastQ_Screen/FastQ-Screen-0.15.3/fastq_screen" #"/path/fastq_screen"
FASTQ_SCREEN_CONF="/storage/shihongjunLab/liulifeng/workflow/soft/FastQ_Screen/FastQ-Screen-0.15.3/fastq_screen.mus.conf" #"/path/fastq_screen.mus.conf"
FASTP_BIN="/storage/shihongjunLab/liulifeng/workflow/soft/bin/fastp" #"/path/fastp"
BOWTIE2_BIN="/soft/bio/bowtie-2.4.2/bowtie2" #"/path/bowtie2"
BOWTIE2_INDEX="/storage/shihongjunLab/liulifeng/database/Genomes/Mus_musculus/Ensembl/GRCm38/Sequence_index/bowtie2/Mus_musculus" #"/path/Mus_musculus"
SAMTOOLS_BIN="/soft/bio/samtools-1.14/bin/samtools" #"/path/samtools"
SAMBAMBA_BIN="/soft/bio/sambamba-0.8.0/sambamba" #"/path/sambamba"
MACS2_BIN="/home/vitoLab/xiongzy12295001/anaconda3/bin/macs2" #"/path/macs2"
BAMCOVERAGE_BIN="/home/vitoLab/xiongzy12295001/anaconda3/bin/bamCoverage" #"/path/bamCoverage"
BLACKLIST="/path/mm10.blacklist.bed" #"/path/mm10.blacklist.bed" 

THREADS=16

# ------------------------- sample list -------------------------
SAMPLES=(ATAC)           # 改成你的样本名

# ============================================================================
for sample in "${SAMPLES[@]}"; do
    echo "[$(date +%F\ %T)] ==================== $sample ===================="

    base="${WORK_DIR}/${sample}"
    r1="${RAW_DIR}/${sample}/${sample}_1.fq.gz"
    r2="${RAW_DIR}/${sample}/${sample}_2.fq.gz"

    if [[ ! -f "$r1" || ! -f "$r2" ]]; then
        echo "ERROR: missing input files for $sample: $r1 / $r2" >&2
        exit 1
    fi

    # ---------------------- 目录 ----------------------
    mkdir -p "$base/01_reads_processing/fastqc" \
             "$base/01_reads_processing/fastq_screen" \
             "$base/01_reads_processing/fastp" \
             "$base/02_mapping_bam_processing/fastq_to_sam" \
             "$base/02_mapping_bam_processing/sam_to_bam" \
             "$base/02_mapping_bam_processing/bam_dup_process" \
             "$base/03_peak_calling/file_bw_bed"

    # ---------------------- 1. QC ----------------------
    "$FASTQC_BIN" -o "$base/01_reads_processing/fastqc" -t "$THREADS" \
        "$r1" "$r2" \
        &> "$base/01_reads_processing/fastqc/fastqc.log"



    "$FASTQ_SCREEN_BIN" --force --Aligner bowtie2 --conf "$FASTQ_SCREEN_CONF" \
        --threads "$THREADS" --subset 1000000 \
        "$r1" "$r2" \
        --outdir "$base/01_reads_processing/fastq_screen" \
        2> "$base/01_reads_processing/fastq_screen/fastq_screen.log"

    # ---------------------- 2. fastp ----------------------
    "$FASTP_BIN" --thread "$THREADS" -z 4 \
        --in1 "$r1" --in2 "$r2" \
        --trim_front1 2 --trim_tail1 0 \
        --trim_front2 2 --trim_tail2 0 \
        --qualified_quality_phred 20 --unqualified_percent_limit 50 \
        --cut_right --cut_window_size 4 --cut_mean_quality 20 \
        --low_complexity_filter \
        --trim_poly_x --trim_poly_g \
        --overrepresentation_analysis \
        --detect_adapter_for_pe \
        --correction \
        --length_required 20 \
        --out1 "$base/01_reads_processing/fastp/${sample}_R1.fastq.gz" \
        --out2 "$base/01_reads_processing/fastp/${sample}_R2.fastq.gz" \
        -j "$base/01_reads_processing/fastp/${sample}.stat.json" \
        -h "$base/01_reads_processing/fastp/${sample}.stat.html" \
        2> "$base/01_reads_processing/fastp/fastp.log"

    # ---------------------- 3. bowtie2 ----------------------
    "$BOWTIE2_BIN" -p "$THREADS" -x "$BOWTIE2_INDEX" \
        -1 "$base/01_reads_processing/fastp/${sample}_R1.fastq.gz" \
        -2 "$base/01_reads_processing/fastp/${sample}_R2.fastq.gz" \
        -S "$base/02_mapping_bam_processing/fastq_to_sam/${sample}.sam" \
        &> "$base/02_mapping_bam_processing/fastq_to_sam/fastq_to_sam.log"

    # ---------------------- 4. sam -> sorted bam ----------------------
    "$SAMTOOLS_BIN" view -b -f 2 -q 30 \
        -o "$base/02_mapping_bam_processing/sam_to_bam/${sample}.unsort_bam" \
        "$base/02_mapping_bam_processing/fastq_to_sam/${sample}.sam"

    "$SAMTOOLS_BIN" sort \
        -o "$base/02_mapping_bam_processing/sam_to_bam/${sample}.sort.bam" \
        "$base/02_mapping_bam_processing/sam_to_bam/${sample}.unsort_bam"

    rm -f "$base/02_mapping_bam_processing/fastq_to_sam/${sample}.sam" \
          "$base/02_mapping_bam_processing/sam_to_bam/${sample}.unsort_bam"

    # ---------------------- 5. MarkDuplicates + 去 MT/Pt + index ----------------------
    "$SAMBAMBA_BIN" markdup --overflow-list-size 600000 \
        --tmpdir="$base/02_mapping_bam_processing" -t "$THREADS" -r \
        "$base/02_mapping_bam_processing/sam_to_bam/${sample}.sort.bam" \
        "$base/02_mapping_bam_processing/bam_dup_process/${sample}.sort.redup.bam" \
        &> "$base/02_mapping_bam_processing/bam_dup_process/pcr_duplicates.log"

    "$SAMTOOLS_BIN" view -h "$base/02_mapping_bam_processing/bam_dup_process/${sample}.sort.redup.bam" \
        | awk 'BEGIN{FS=OFS="\t"} /^@/{print; next} $3!="MT"' \
        | "$SAMTOOLS_BIN" view -bS \
            -o "$base/02_mapping_bam_processing/bam_dup_process/${sample}.final.bam"

    "$SAMTOOLS_BIN" index "$base/02_mapping_bam_processing/bam_dup_process/${sample}.final.bam"

    # ---------------------- 6. MACS2 peak calling ----------------------
    "$MACS2_BIN" callpeak \
        -t "$base/02_mapping_bam_processing/bam_dup_process/${sample}.final.bam" \
        -f BAMPE -g mm -n "$sample" \
        --outdir "$base/03_peak_calling" \
        2> "$base/03_peak_calling/${sample}.macs2.log"



    # ---------------------- 7. blacklist filter ----------------------
    bedtools intersect -v \
        -a "$base/03_peak_calling/${sample}.peaks.narrowPeak" \
        -b "$BLACKLIST" \
        > "$base/03_peak_calling/${sample}.peaks.narrowPeak.tmp"
    mv "$base/03_peak_calling/${sample}.peaks.narrowPeak.tmp" \
       "$base/03_peak_calling/${sample}.peaks.narrowPeak"


    # ---------------------- 8. bigwig + bed ----------------------
    "$BAMCOVERAGE_BIN" --normalizeUsing CPM \
        -b "$base/02_mapping_bam_processing/bam_dup_process/${sample}.final.bam" \
        -o "$base/03_peak_calling/file_bw_bed/${sample}.final.bw" \
        &> "$base/03_peak_calling/file_bw_bed/bw.log"
done

echo "[$(date +%F\ %T)] All samples done."