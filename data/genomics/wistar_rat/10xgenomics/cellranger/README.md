# Cellranger test datasets

This folder contains test datasets from 10X Genomics for reference and testing its proprietary [Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) software.

# Data sources

| Directory | Source | Technologies | Cellranger calls |
| --------- | ------ | ------------ | ---------------- |
| 10k_pbmc_ocm | [10k Wistar Rat PBMCs Multiplex Sample](https://www.10xgenomics.com/datasets/10k_Wistar_Rat_PBMCs_Multiplex_3p_gem-x_Universal_OCM) | GEX, Cell Multiplexing | `multi` |

# Subsampling

The original datasets contain FASTQ files that are too large to store here.
Unless stated otherwise, FASTQs were naively subsampled to 10,000 reads by reading the first 40,000 lines of each FASTQ file (4 lines per read).

# Directory structure

```bash
.
├── 10k_pbmc_ocm
│   ├── 10k_pbmc_cmo_config.csv
│   ├── fastqs
│   │   └── gex
│   │       ├── 10k_Wistar_Rat_PBMCs_Multiplex_3p_gem-x_Universal_OCM_S1_L001_R1_001.fastq.gz
│   │       ├── 10k_Wistar_Rat_PBMCs_Multiplex_3p_gem-x_Universal_OCM_S1_L001_R2_001.fastq.gz
│   │       ├── 10k_Wistar_Rat_PBMCs_Multiplex_3p_gem-x_Universal_OCM_S1_L002_R1_001.fastq.gz
│   │       ├── 10k_Wistar_Rat_PBMCs_Multiplex_3p_gem-x_Universal_OCM_S1_L002_R2_001.fastq.gz
│   │       ├── 10k_Wistar_Rat_PBMCs_Multiplex_3p_gem-x_Universal_OCM_S1_L003_R1_001.fastq.gz
│   │       ├── 10k_Wistar_Rat_PBMCs_Multiplex_3p_gem-x_Universal_OCM_S1_L003_R2_001.fastq.gz
│   │       ├── 10k_Wistar_Rat_PBMCs_Multiplex_3p_gem-x_Universal_OCM_S1_L004_R1_001.fastq.gz
│   │       ├── 10k_Wistar_Rat_PBMCs_Multiplex_3p_gem-x_Universal_OCM_S1_L004_R2_001.fastq.gz
│   │       ├── 10k_Wistar_Rat_PBMCs_Multiplex_3p_gem-x_Universal_OCM_S1_L005_R1_001.fastq.gz
│   │       ├── 10k_Wistar_Rat_PBMCs_Multiplex_3p_gem-x_Universal_OCM_S1_L005_R2_001.fastq.gz
│   │       ├── 10k_Wistar_Rat_PBMCs_Multiplex_3p_gem-x_Universal_OCM_S1_L006_R1_001.fastq.gz
│   │       └── 10k_Wistar_Rat_PBMCs_Multiplex_3p_gem-x_Universal_OCM_S1_L006_R2_001.fastq.gz
│   └── README.md
└── README.md
```
