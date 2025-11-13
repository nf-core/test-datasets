# Rare Disease Test Datasets  

This repository contains subsampled long-read sequencing datasets tailored for rare disease analysis.  


---

## Contents  

- `bam_pass/` – subsampled aligned unmapped BAM files for variant calling tests from HG002 basecalled data
- `fastq_file/` - subsamples fastq file from HG002 basecalled data
- `spectre/` – VCF files and BED regions for whole-genome CNV testing  
- `straglr/` – Chromosome 22 STR test regions  
- `test.exclude.bed` – CNV test exclude regions  
- `reference/` – reduced human genome references   
- `samplesheet_*.csv` – metadata for pipeline test runs  

---

## Sample Overview  

| Sample ID | File type   | Size (approx.) | Purpose                                      |  
|-----------|-------------|----------------|----------------------------------------------|  
| Test      | BAM         | ~100 MB        | End-to-end pipeline testing from alignment (minimap2) through variant analysis | 
| Reference | FASTA / BED | <5 MB          | Subset references (Chromosome 22) for rare disease test runs                   |  

---

## Usage  

These datasets are intended for automated testing of long-read rare disease pipeline (https://github.com/nf-core/longraredisease).  

The data in this repository will be used to test the pipeline starting from unaligned BAM files (using minimap2).  
The associated parameters and settings to run the pipeline can be found in the **test.config** file.  

Example run:  

```bash
nextflow run nf-core/nanoraredx -profile test,docker
