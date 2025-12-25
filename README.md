# Small RNA‑seq (miRNA‑seq) pipeline for ulcerative colitis blood

## Overview

This repository contains an end‑to‑end **small RNA‑seq (miRNA‑seq)** analysis pipeline for human peripheral blood samples from ulcerative colitis (UC) patients and healthy controls. It uses a small subset of the GEO dataset **GSE169570** (SRA project **PRJNA717025**) and standard tools (cutadapt, bowtie, samtools, DESeq2). The pipeline starts from raw FASTQ files and produces miRNA‑level counts and differential expression results. The main goal is to demonstrate a clean, reproducible NGS workflow on real data.

## Dataset

**Study:** GSE169570 – small RNA‑seq of peripheral blood  
**Subset used here (3 UC, 1 control):**

| Sample ID   | Condition          |
|------------|--------------------|
| SRR14062898 | Ulcerative colitis |
| SRR14062899 | Ulcerative colitis |
| SRR14062940 | Ulcerative colitis |
| SRR14062947 | Healthy control    |

Raw data are downloaded from SRA (not stored in this repo). Download commands are in `scripts/01_download_data.sh`.

## Pipeline

1. Download raw small RNA‑seq reads from SRA for the four samples.  
2. Trim Illumina TruSeq Small RNA 3′ adapter and keep 15–30 nt inserts using **cutadapt**.  
3. Download `mature.fa` from **miRBase**, extract human (`hsa-`) mature miRNAs, and build a **bowtie** index (`hsa_mature`).  
4. Align trimmed reads to `hsa_mature` with bowtie, then convert SAM → sorted BAM and index using **samtools**.  
5. Generate per‑miRNA counts with `samtools idxstats` and combine them into `results/miRNA_counts.tsv`.  
6. Run **DESeq2** in R with design `~ condition` (UC vs control) to obtain `results/DESeq2_miRNA_results.csv`.

## How to run

**Requirements**

- Linux + bash  
- Tools: `cutadapt`, `bowtie` (v1), `samtools`  
- R (≥ 4.0) with the `DESeq2` package

**Commands**

From the repo root:

bash scripts/01_download_data.sh
bash scripts/02_trim_cutadapt.sh
bash scripts/03_build_miRBase_index.sh
bash scripts/04_align_bowtie.sh
bash scripts/05_counts_idxstats.sh
Rscript scripts/06_deseq2_analysis.R

text

Key outputs:

- `results/miRNA_counts.tsv` – miRNA × sample count matrix  
- `results/DESeq2_miRNA_results.csv` – differential expression results (UC vs control)

## Results & limitations

On this small subset (3 UC vs 1 control), DESeq2 finds a few miRNAs with notable log2 fold changes between UC and control, but only one miRNA reaches an adjusted p‑value around 0.1. Because there is only one control sample, these results are **exploratory** and not suitable for strong biological conclusions. The primary purpose of this project is to demonstrate a reproducible small RNA‑seq / miRNA‑seq pipeline, including NGS preprocessing and DE analysis, on real human data.


## Results (short summary)

**Alignment / QC**

- After adapter trimming and length filtering, between ~3.8–5.5 million reads per sample remained.  
- Bowtie alignment to human mature miRNAs (miRBase) achieved **~70–81%** mapped reads across the four samples.  
- This mapping rate is consistent with expectations for peripheral blood small RNA‑seq libraries.

**miRNA counts and DESeq2**

- `samtools idxstats` produced non‑zero counts for **407** miRNAs in total.  
- DESeq2 was run with design `~ condition` (UC vs control) after filtering out very low‑count miRNAs.  
- On this subset (3 UC vs 1 control), DESeq2 reported:
  - **1 miRNA** with adjusted p‑value < 0.1 (up‑regulated in UC).  
  - Example: `hsa-miR-3940-3p` with log2 fold change ≈ **+3.6** (higher in UC) and padj ≈ **0.06**.  
  - Several other miRNAs show notable log2 fold changes but have padj ≥ 0.1 due to the very small sample size.

**Interpretation**

- Alignment and counting confirm that the pipeline correctly captures known blood miRNAs with reasonable mapping rates.  
- The differential expression results are **exploratory only**, because having 3 UC vs 1 control provides limited power and unstable variance estimates for the control group.  
- The main value of this project is demonstrating a reproducible small RNA‑seq / miRNA‑seq pipeline (download → trim → align → count → DESeq2), not delivering final biomarker claims.




