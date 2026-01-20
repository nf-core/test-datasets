# Cellranger test datasets

This folder contains test datasets from 10X Genomics for reference and testing its proprietary [Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) software.

# Data sources

| Directory | Source | Technologies | Cellranger calls |
| --------- | ------ | ------------ | ---------------- |
|5k_cmvpos_tcells | [Integrated GEX, TotalSeq™-C, and TCR Analysis of Chromium Connect Generated Library from 5k CMV+ T cells](https://www.10xgenomics.com/resources/datasets/integrated-gex-totalseqc-and-tcr-analysis-of-connect-generated-library-from-5k-cmv-t-cells-2-standard) | GEX, TCR, Antibody Capture | `count`, `vdj`, `multi` | 
|10k_pbmc | [Human PBMC from a Healthy Donor, 10k cells - multi (v2)](https://www.10xgenomics.com/resources/datasets/human-pbmc-from-a-healthy-donor-10-k-cells-multi-v-2-2-standard-5-0-0) | GEX, Fixed RNA Profiling, V(D)J-B, V(D)J-T, Antibody Capture | `count`, `vdj`, `multi` |
| 10k_pbmc_cmo | [10k Human PBMCs Stained with TotalSeq™-B Human Universal Cocktail, Singleplex Sample](https://www.10xgenomics.com/resources/datasets/10k-human-pbmcs-stained-with-totalseq-b-human-universal-cocktail-singleplex-sample-1-standard) | GEX, Cell Multiplexing | `count`, `multi` |
| 4plex_scFFPE | [Mixture of Healthy and Cancer FFPE Tissues Dissociated using Miltenyi FFPE Tissue Dissociation Kit, Multiplexed Samples, 4 Probe Barcodes](https://www.10xgenomics.com/datasets/mixture-of-healthy-and-cancer-ffpe-tissues-dissociated-using-miltenyi-ffpe-tissue-dissociation-kit-multiplexed-samples-4-probe-barcodes-1-standard) | GEX, FFPE, Cell Multiplexing | `multi` |
| sc3_v3_5k_a549_gex_crispr | [5k A549, Lung Carcinoma Cells, No Treatment Transduced with a CRISPR Pool](https://www.10xgenomics.com/datasets/5-k-a-549-lung-carcinoma-cells-no-treatment-transduced-with-a-crispr-pool-3-1-standard-6-0-0) | GEX, CRISPR | `count`, `multi` |
| sc5p_v2_hs_PBMC_1k_bcr | [AIRR rearrangement (TSV) from Single Cell Immune Profiling Dataset](https://support.10xgenomics.com/single-cell-vdj/datasets/4.0.0/sc5p_v2_hs_PBMC_1k) | BCR, AIRR rearrangement | `vdj` |
| singleplex_flex | [Human Kidney Nuclei GEM-X Flex Gene Expression](https://www.10xgenomics.com/datasets/Human_Kidney_4k_GEM-X_Flex) | GEX, Flex, Singleplex | `multi` |

# Subsampling

The original datasets contain FASTQ files that are too large to store here.
Unless stated otherwise, FASTQs were naively subsampled to 10,000 reads by reading the first 40,000 lines of each FASTQ file (4 lines per read).

# Directory structure

```bash
.
├── 10k_pbmc
│   ├── 10k_pbmc_config.csv
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
│   ├── README.md
│   └── sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t_feature_ref.csv
├── 10k_pbmc_cmo
│   ├── 10k_pbmc_cmo_config.csv
│   ├── 10k_pbmc_cmo_count_feature_reference.csv
│   ├── fastqs
│   │   ├── cmo
│   │   │   ├── subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_multiplexing_capture_S1_L001_R1_001.fastq.gz
│   │   │   └── subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_multiplexing_capture_S1_L001_R2_001.fastq.gz
│   │   ├── gex_1
│   │   │   ├── subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_gex_S2_L001_R1_001.fastq.gz
│   │   │   └── subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_gex_S2_L001_R2_001.fastq.gz
│   │   └── gex_2
│   │       ├── subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_2_gex_S1_L001_R1_001.fastq.gz
│   │       └── subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_2_gex_S1_L001_R2_001.fastq.gz
│   └── README.md
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
│   ├── fastqs
│   │   ├── ab
│   │   │   ├── subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R1_001.fastq.gz
│   │   │   └── subsampled_5k_human_antiCMV_T_TBNK_connect_AB_S2_L004_R2_001.fastq.gz
│   │   ├── gex_1
│   │   │   ├── subsampled_5k_human_antiCMV_T_TBNK_connect_GEX_1_S1_L001_R1_001.fastq.gz
│   │   │   └── subsampled_5k_human_antiCMV_T_TBNK_connect_GEX_1_S1_L001_R2_001.fastq.gz
│   │   └── vdj
│   │       ├── subsampled_5k_human_antiCMV_T_TBNK_connect_VDJ_S1_L001_R1_001.fastq.gz
│   │       └── subsampled_5k_human_antiCMV_T_TBNK_connect_VDJ_S1_L001_R2_001.fastq.gz
│   └── README.md
├── hashing_demultiplexing
│   ├── 438-21-raw_feature_bc_matrix.h5
│   ├── 438_21_raw_HTO.csv
│   ├── hashsolo_anndata.h5ad
│   ├── hto
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   ├── htodemux.rds
│   ├── hto.tar.gz
│   ├── README.md
│   ├── rna
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   └── rna.tar.gz
├── README.md
├── references
│   ├── flex
│   │   └── Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.chr22.csv
│   ├── README.md
│   └── vdj
│       └── refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0
│           ├── fasta
│           │   ├── regions.fa
│           │   └── supp_regions.fa
│           └── reference.json
├── sc3_v3_5k_a549_gex_crispr
│   ├── fastqs
│   │   ├── crispr
│   │   │   ├── subsampled_SC3_v3_NextGem_DI_CRISPR_A549_5K_crispr_S4_L001_R1_001.fastq.gz
│   │   │   ├── subsampled_SC3_v3_NextGem_DI_CRISPR_A549_5K_crispr_S4_L001_R2_001.fastq.gz
│   │   │   ├── subsampled_SC3_v3_NextGem_DI_CRISPR_A549_5K_crispr_S4_L002_R1_001.fastq.gz
│   │   │   └── subsampled_SC3_v3_NextGem_DI_CRISPR_A549_5K_crispr_S4_L002_R2_001.fastq.gz
│   │   └── gex
│   │       ├── subsampled_SC3_v3_NextGem_DI_CRISPR_A549_5K_gex_S5_L001_R1_001.fastq.gz
│   │       ├── subsampled_SC3_v3_NextGem_DI_CRISPR_A549_5K_gex_S5_L001_R2_001.fastq.gz
│   │       ├── subsampled_SC3_v3_NextGem_DI_CRISPR_A549_5K_gex_S5_L002_R1_001.fastq.gz
│   │       └── subsampled_SC3_v3_NextGem_DI_CRISPR_A549_5K_gex_S5_L002_R2_001.fastq.gz
│   ├── README.md
│   ├── reference
│   │   ├── genes_chr1_32292083_32686211.gtf
│   │   └── genome_chr1_32292083_32686211.fa
│   ├── SC3_v3_NextGem_DI_CRISPR_A549_5K_Multiplex_config.csv
│   ├── SC3_v3_NextGem_DI_CRISPR_A549_5K_Multiplex_count_feature_reference_chr1_32292083_32686211.csv
│   └── SC3_v3_NextGem_DI_CRISPR_A549_5K_Multiplex_count_feature_reference.csv
├── sc5p_v2_hs_PBMC_1k_bcr
│   └── sc5p_v2_hs_PBMC_1k_b_airr_rearrangement.tsv
└── singleplex_flex
    ├── Human_Kidney_GEM-X_Flex_S1_L001_R1_001.subsampled.fastq.gz
    └── Human_Kidney_GEM-X_Flex_S1_L001_R2_001.subsampled.fastq.gz
```
