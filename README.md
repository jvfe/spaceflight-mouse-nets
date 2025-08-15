# README

## 1. Data sources

This task is based on publicly available sequencing data from a study of the contrasting patterns of gene interactions and radiation dose-dependent effects in multiple tissues of spaceflight mice. The dataset includes samples under flight and ground control conditions and was originally sequenced using Illumina HiSeq.

---

## 2. How to download

The data can be downloaded from the NASA Open Science Data Repository (OSDR) under accession number OSD-47 (samples FLT1 and GC2).

---

## 3. Subsampling

To create a focused dataset for this workflow, the original data was subsampled to only include reads that align to chromosomes 17, 18, and 19 of the mouse GRCm39 reference genome.

- Align original FASTQ files to a reduced GRCm39 reference containing only data for chr17, chr18, and chr19.

---

## 4. How the workflow works

The workflow files are stored in workflow/. All scripts should be executed from the spaceflight-mouse directory.

---

### Step 1 – QC, Align and Quantify

**Purpose:** Align subsampled FASTQ reads to the reference (GRCm39 chr17-19) and quantify gene expression.
**Script:** workflow/run_pipeline.sh
**Inputs:** FASTQ files from data/raw_data
**Outputs:** Alignment files (BAM) and raw gene counts for each sample.
**Command:**

```bash
bash workflow/run_pipeline.sh
```

---

### Step 2 Gather Counts

**Purpose:** Combine the individual count files into a single count matrix.
**Script:** workflow/gather_counts.R
**Inputs:** Raw count files generated from Step 1.
**Outputs:** A single, unified count matrix file.
**Command:**

```bash
Rscript workflow/gather_counts.R
```

---

### Step 3 – Generate Single-Sample Networks and Enrichment

**Purpose:** Create a gene co-expression network for each individual sample and perform functional enrichment analysis.
**Script:** workflow/lioness.R
**Inputs:** The unified count matrix from Step 2 as well as the provided metadata table.
**Outputs:** Significant differentially interacting genes (DIGs) and enrichment analysis results.
**Command:**

```bash
Rscript workflow/lioness.R
```
