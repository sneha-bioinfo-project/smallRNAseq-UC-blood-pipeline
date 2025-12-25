#!/usr/bin/env bash
set -euo pipefail

########## CONFIG ##########
RAW_DIR="../raw_fastq"
TRIM_DIR="../trimmed_fastq"
FASTQC_RAW_DIR="../fastqc/raw"
REF_DIR="../ref"
ALIGN_DIR="../align_bowtie"
BAM_DIR="../align_bowtie_bam"
REPORT_DIR="../reports/cutadapt"
RESULTS_DIR="../results"
PLOTS_DIR="../results/plots"

# Typical small RNA 3' adapter; change if your kit uses a different one. [web:608][web:576]
ADAPTER="TGGAATTCTCGGGTGCCAAGG"
MINLEN=15   # keep reads >= 15 nt for miRNA. [web:606]
MIRBASE_FA="${REF_DIR}/hsa_mature.fa"   # mature miRNA fasta from miRBase
INDEX_PREFIX="${REF_DIR}/hsa_mature"    # bowtie index prefix

mkdir -p "${TRIM_DIR}" "${FASTQC_RAW_DIR}" "${ALIGN_DIR}" "${BAM_DIR}" \
         "${REPORT_DIR}" "${RESULTS_DIR}" "${PLOTS_DIR}"

########## STEP 1: FastQC on raw reads ##########
echo "=== Step 1: FastQC on raw FASTQ ==="
fastqc -t 4 -o "${FASTQC_RAW_DIR}" "${RAW_DIR}"/*.fastq.gz   # standard QC step. [web:608][web:602]

########## STEP 2: Adapter trimming with cutadapt ##########
echo "=== Step 2: Adapter trimming with cutadapt ==="
for fq in "${RAW_DIR}"/*.fastq.gz
do
    base=$(basename "${fq}" .fastq.gz)
    cutadapt \
        -a "${ADAPTER}" \
        -m "${MINLEN}" \
        -o "${TRIM_DIR}/${base}_trimmed.fastq.gz" \
        "${fq}" \
        > "${REPORT_DIR}/${base}_cutadapt.txt"
done
# cutadapt is widely used for small RNA adapter removal. [web:603][web:608]

########## STEP 3: Build Bowtie index (miRBase) ##########
echo "=== Step 3: Build Bowtie index for miRBase (if missing) ==="
if [ ! -f "${INDEX_PREFIX}.1.ebwt" ]; then
    bowtie-build "${MIRBASE_FA}" "${INDEX_PREFIX}"   # Bowtie1 is often used for miRNA. [web:612][web:609]
fi

########## STEP 4: Align trimmed reads with Bowtie ##########
echo "=== Step 4: Align trimmed reads to miRBase ==="
for fq in "${TRIM_DIR}"/*_trimmed.fastq.gz
do
    base=$(basename "${fq}" _trimmed.fastq.gz)

    # Map short miRNA reads with Bowtie v1. [web:576][web:612]
    bowtie \
        -v 1 -m 10 --best --strata \
        -q "${INDEX_PREFIX}" \
        -U "${fq}" \
        -S "${ALIGN_DIR}/${base}.sam"

    samtools view -bS "${ALIGN_DIR}/${base}.sam" \
    | samtools sort -o "${BAM_DIR}/${base}.sorted.bam" -

    samtools index "${BAM_DIR}/${base}.sorted.bam"
done

########## STEP 5: Build miRNA counts matrix (idxstats) ##########
echo "=== Step 5: Build counts matrix from BAM ==="
OUT="${RESULTS_DIR}/miRNA_counts.tsv"

samples=()
for bam in "${BAM_DIR}"/*.sorted.bam
do
    samples+=("$(basename "${bam}" .sorted.bam)")
done

{
    printf "miRNA"
    for s in "${samples[@]}"; do printf "\t%s" "${s}"; done
    printf "\n"

    first_bam="${BAM_DIR}/${samples[0]}.sorted.bam"
    cut -f1 <(samtools idxstats "${first_bam}") \
        | while read -r mir; do
            [ "${mir}" = "*" ] && continue
            printf "%s" "${mir}"
            for s in "${samples[@]}"; do
                count=$(samtools idxstats "${BAM_DIR}/${s}.sorted.bam" \
                        | awk -v id="${mir}" '$1==id {print $3}')
                printf "\t%s" "${count:-0}"
            done
            printf "\n"
        done
} > "${OUT}"
# Expression matrices for small RNA are then used for differential analysis. [web:616]

echo "=== Finished alignment + counts. Now run DESeq2 R script ==="

