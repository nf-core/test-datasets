# 5k_cmvpos_tcells

Human T cells were isolated from a patient who was previously diagnosed with cytomegalovirus (lot 5147JN21).
Gene Expression, TCR, and Antibody Capture libraries were generated from a single Chromium Connect channel.
Antibody libraries were generated using a combination of a BioLegend TotalSeq™-C TBNK panel, a CMV dextramer, and a non-binding control dextramer.
The (original) TCR library was downsampled to 5,000 reads per cell.
Libraries were prepared following the Chromium Next GEM Automated Single Cell 5' Reagent Kits v2 User Guide (CG0000507).

Data were sequenced on Illumina NovaSeq with 28 bp read 1, 90 bp read 2, 10 bp i5 sample barcode, and 10 bp i7 sample barcode.

## Download original input data

```bash
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-vdj/6.1.2/5k_human_antiCMV_T_TBNK_connect_Multiplex/5k_human_antiCMV_T_TBNK_connect_Multiplex_fastqs.tar
wget https://cf.10xgenomics.com/samples/cell-vdj/6.1.2/5k_human_antiCMV_T_TBNK_connect_Multiplex/5k_human_antiCMV_T_TBNK_connect_Multiplex_config.csv
wget https://cf.10xgenomics.com/samples/cell-vdj/6.1.2/5k_human_antiCMV_T_TBNK_connect_Multiplex/5k_human_antiCMV_T_TBNK_connect_Multiplex_count_feature_reference.csv
```

## Original `cellranger multi` config

```bash
[gene-expression]
ref,/path/to/references/refdata-gex-GRCh38-2020-A
expect-cells,5000

[vdj]
ref,/path/to/vdj_references/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0

[feature]
ref,/path/to/feature_references/BioL_TotalseqC_Immudex_CMV_TCR_1_Dextramer.csv

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
5k_human_antiCMV_T_TBNK_connect_GEX_1,/path/to/fastqs/5k_human_antiCMV_T_TBNK_connect/gex_1,1-4,gex,gene expression,
5k_human_antiCMV_T_TBNK_connect_GEX_2,/path/to/fastqs/5k_human_antiCMV_T_TBNK_connect/gex_2,3-4,gex,gene expression,
5k_human_antiCMV_T_TBNK_connect_AB,/path/to/fastqs/5k_human_antiCMV_T_TBNK_connect/ab,4,ab,antibody capture,
5k_human_antiCMV_T_TBNK_connect_VDJ,/path/to/fastqs/5k_human_antiCMV_T_TBNK_connect/vdj,1,vdj,vdj-t,0.422012153950034
```
