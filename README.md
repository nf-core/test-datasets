# 🧬 Rare Disease Test Datasets

This repository provides **subsampled Oxford Nanopore long-read sequencing datasets** derived from **HG002**, designed for **testing and validation** of long-read rare disease analysis pipelines such as [nf-core/longraredisease](https://github.com/nf-core/longraredisease).

All datasets are restricted to **chromosome 22 (first 50 Mb)** to minimise file sizes and speed up automated test runs.

---

## 📂 Repository Contents

| Folder / File       | Description                                                                                |
| ------------------- | ------------------------------------------------------------------------------------------ |
| `ubam_file/`        | Subsampled **unmapped BAM** files (uBAMs) for testing variant calling from unaligned data. |
| `fastq_file/`       | Subsampled **FASTQ** file generated from HG002 basecalled reads.                           |
| `spectre/`          | Example **VCF** and **BED** files for CNV detection testing with _Spectre_.                |
| `straglr/`          | **STR test regions** (chromosome 22) for _STRaglr_ validation.                             |
| `hificnv/`          | **Exclude BED** regions used for chromosome 22 CNV benchmarking.                           |
| `reference/`        | Reduced **human genome reference**, containing only chromosome 22 (GRCh38).                |
| `samplesheet_*.csv` | Example **sample metadata** for automated pipeline test runs.                              |

---

## 🧪 Sample Overview

| Column                                    | Description                           |
| ----------------------------------------- | ------------------------------------- |
| `sample_id`                               | Unique identifier for the test sample |
| `input_type`                              | Input data type (FASTQ, BAM, etc.)    |
| `file_path`                               | Direct download link to test data     |
| `hpo_terms`                               | Associated HPO phenotype terms        |
| `sex`                                     | Biological sex                        |
| `family_id`, `maternal_id`, `paternal_id` | Family metadata                       |

Example entry:

```
sample_id,input_type,file_path,hpo_terms,sex,family_id,maternal_id,paternal_id
test,fastq,https://raw.githubusercontent.com/nourmahfel/test-datasets/longraredisease/fastq_file/hg002_subset.fastq.gz,HP:0002721;HP:0002110;HP:0500093;HP:0000717;HP:0001263;HP:0001763;HP:0003298;HP:0002857;HP:0001382,F,family_21,null,null
```

---

## ⚙️ Usage

These datasets are intended for **automated pipeline testing**, enabling quick validation of the full _long-read rare disease analysis_ workflow — from unaligned reads through to variant calling and annotation.

Example Nextflow test run:

```bash
nextflow run nf-core/nanoraredx -profile test,docker
```

The repository includes a `test.config` file containing preset paths and parameters used for CI and development validation.

---

## 🧩 Data Generation Workflow

The following steps describe how each dataset was created from **HG002 Oxford Nanopore sequencing data**.

### 1️⃣ Extract 50 Mb region from chromosome 22

A compact subset was created to minimise storage and runtime while preserving data realism:

```bash
samtools view -b calls.sorted.bam chr22:1-50000000 > chr22_50mb.bam
samtools index chr22_50mb.bam
```

This produced a **50 Mb** region representing chromosome 22 (`chr22_50mb.bam`).

---

### 2️⃣ Generate FASTQ file

The BAM file was converted to FASTQ format to simulate basecalled reads:

```bash
samtools fastq chr22_50mb.bam > hg002_subset.fastq
gzip hg002_subset.fastq
```

Output:

- `hg002_subset.fastq.gz` → FASTQ dataset for testing pipeline entry from raw reads.

---

### 3️⃣ Create unmapped BAM (uBAM)

To test the alignment and variant calling stages from unaligned data, an **unmapped BAM** version was generated:

```bash
samtools view -h chr22_50mb.bam | awk '$3=="*" || /^@/' | samtools view -b -o hg002_subset.ubam
```

This file retains read names, qualities, and tags but removes alignment fields (RNAME, POS, CIGAR, etc.).

Output:

- `hg002_subset.ubam` → unaligned BAM file suitable for pipeline tests starting from mapping.

---

## 📦 Summary of Derived Outputs

| File                             | Description                        | Source           |
| -------------------------------- | ---------------------------------- | ---------------- |
| `hg002_subset.fastq.gz`          | Subsampled FASTQ (50 Mb region)    | `chr22_50mb.bam` |
| `hg002_subset.ubam`              | Unmapped BAM for alignment testing | `chr22_50mb.bam` |
| `spectre/*.vcf`, `spectre/*.bed` | CNV test data                      | chr22 subset     |
| `straglr/*.bed`                  | STR test regions                   | chr22 subset     |
| `reference/chr22.fasta`          | Reduced genome reference           | GRCh38           |

---

## 📄 License and Attribution

Data derived from **HG002** (Genome in a Bottle Consortium).  
Please cite **GIAB** and relevant tools when reusing or redistributing these datasets.
