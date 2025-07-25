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
├── 10k_pbmc
│   ├── 10k_pbmc_config.csv
│   ├── README.md
│   ├── fastqs
│   │   ├── 5gex
│   │   │   ├── 5fb
│   │   │   │   ├── subsampled_sc5p_v2_hs_PBMC_10k_5fb_S1_L001_R1_001.fastq.gz
│   │   │   │   └── subsampled_sc5p_v2_hs_PBMC_10k_5fb_S1_L001_R2_001.fastq.gz
│   │   │   └── 5gex
│   │   │       ├── subsampled_sc5p_v2_hs_PBMC_10k_5gex_S1_L001_R1_001.fastq.gz
│   │   │       └── subsampled_sc5p_v2_hs_PBMC_10k_5gex_S1_L001_R2_001.fastq.gz
│   │   ├── bcell
│   │   │   ├── subsampled_sc5p_v2_hs_PBMC_10k_b_S1_L001_R1_001.fastq.gz
│   │   │   └── subsampled_sc5p_v2_hs_PBMC_10k_b_S1_L001_R2_001.fastq.gz
│   │   └── tcell
│   │       ├── subsampled_sc5p_v2_hs_PBMC_10k_t_S1_L001_R1_001.fastq.gz
│   │       └── subsampled_sc5p_v2_hs_PBMC_10k_t_S1_L001_R2_001.fastq.gz
│   └── sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t_feature_ref.csv
├── 10k_pbmc_cmo
│   ├── 10k_pbmc_cmo_config.csv
│   ├── 10k_pbmc_cmo_count_feature_reference.csv
│   ├── README.md
│   └── fastqs
│       ├── cmo
│       │   ├── subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_multiplexing_capture_S1_L001_R1_001.fastq.gz
│       │   └── subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_multiplexing_capture_S1_L001_R2_001.fastq.gz
│       ├── gex_1
│       │   ├── subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_gex_S2_L001_R1_001.fastq.gz
│       │   └── subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_gex_S2_L001_R2_001.fastq.gz
│       └── gex_2
│           ├── subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_2_gex_S1_L001_R1_001.fastq.gz
│           └── subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_2_gex_S1_L001_R2_001.fastq.gz
├── 4plex_scFFPE
│   ├── 4plex_human_liver_colorectal_ovarian_panc_scFFPE_multiplex_S1_L001_R1_001.subsampled.fastq.gz
│   ├── 4plex_human_liver_colorectal_ovarian_panc_scFFPE_multiplex_S1_L001_R2_001.subsampled.fastq.gz
│   ├── 4plex_human_liver_colorectal_ovarian_panc_scFFPE_multiplex_S1_L002_R1_001.subsampled.fastq.gz
│   ├── 4plex_human_liver_colorectal_ovarian_panc_scFFPE_multiplex_S1_L002_R2_001.subsampled.fastq.gz
│   ├── 4plex_human_liver_colorectal_ovarian_panc_scFFPE_multiplex_S1_L003_R1_001.subsampled.fastq.gz
│   ├── 4plex_human_liver_colorectal_ovarian_panc_scFFPE_multiplex_S1_L003_R2_001.subsampled.fastq.gz
│   ├── 4plex_human_liver_colorectal_ovarian_panc_scFFPE_multiplex_S1_L004_R1_001.subsampled.fastq.gz
│   └── 4plex_human_liver_colorectal_ovarian_panc_scFFPE_multiplex_S1_L004_R2_001.subsampled.fastq.gz
├── 5k_cmvpos_tcells
│   ├── 5k_cmvpos_tcells_config.csv
│   ├── 5k_human_antiCMV_T_TBNK_connect_Multiplex_count_feature_reference.csv
│   ├── README.md
│   └── fastqs
│       ├── ab
│       │   ├── subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R1_001.fastq.gz
│       │   └── subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R2_001.fastq.gz
│       ├── gex_1
│       │   ├── subsampled_5k_human_antiCMV_T_TBNK_connect_GEX_1_S1_L001_R1_001.fastq.gz
│       │   └── subsampled_5k_human_antiCMV_T_TBNK_connect_GEX_1_S1_L001_R2_001.fastq.gz
│       └── vdj
│           ├── subsampled_5k_human_antiCMV_T_TBNK_connect_VDJ_S1_L001_R1_001.fastq.gz
│           └── subsampled_5k_human_antiCMV_T_TBNK_connect_VDJ_S1_L001_R2_001.fastq.gz
├── README.md
├── hashing_demultiplexing
│   ├── 438-21-raw_feature_bc_matrix.h5
│   ├── 438_21_raw_HTO.csv
│   ├── README.md
│   ├── hto
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   ├── hto.tar.gz
│   ├── rna
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   └── rna.tar.gz
├── references
│   ├── README.md
│   └── vdj
│       └── refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0
│           ├── fasta
│           │   ├── regions.fa
│           │   └── supp_regions.fa
│           └── reference.json
├── sc3_v3_5k_a549_gex_crispr
│   ├── fastqs
│   │   ├── crispr
│   │   │   ├── subsampled_SC3_v3_NextGem_DI_CRISPR_A549_5K_crispr_S4_L001_R1_001.fastq.gz
│   │   │   ├── subsampled_SC3_v3_NextGem_DI_CRISPR_A549_5K_crispr_S4_L001_R2_001.fastq.gz
│   │   │   ├── subsampled_SC3_v3_NextGem_DI_CRISPR_A549_5K_crispr_S4_L002_R1_001.fastq.gz
│   │   │   └── subsampled_SC3_v3_NextGem_DI_CRISPR_A549_5K_crispr_S4_L002_R2_001.fastq.gz
│   │   └── gex
│   │       ├── subsampled_SC3_v3_NextGem_DI_CRISPR_A549_5K_gex_S5_L001_R1_001.fastq.gz
│   │       ├── subsampled_SC3_v3_NextGem_DI_CRISPR_A549_5K_gex_S5_L001_R2_001.fastq.gz
│   │       ├── subsampled_SC3_v3_NextGem_DI_CRISPR_A549_5K_gex_S5_L002_R1_001.fastq.gz
│   │       └── subsampled_SC3_v3_NextGem_DI_CRISPR_A549_5K_gex_S5_L002_R2_001.fastq.gz
│   ├── README.md
│   ├── reference
│   │   ├── genes_chr1_32292083_32686211.gtf
│   │   └── genome_chr1_32292083_32686211.fa
│   ├── SC3_v3_NextGem_DI_CRISPR_A549_5K_Multiplex_config.csv
│   ├── SC3_v3_NextGem_DI_CRISPR_A549_5K_Multiplex_count_feature_reference_chr1_32292083_32686211.csv
│   └── SC3_v3_NextGem_DI_CRISPR_A549_5K_Multiplex_count_feature_reference.csv
└── sc5p_v2_hs_PBMC_1k_bcr
    └── sc5p_v2_hs_PBMC_1k_b_airr_rearrangement.tsv

```
