# 10k Human PBMCs Multiplexed, 2 CMOs

Peripheral blood mononuclear cells (PBMCs) from a healthy 19 year old female donor were obtained by 10x Genomics from IQ Biosciences and multiplexed at equal proportions with 2 CMOs.

The config.csv input file was submitted with two sample IDs (correlating to the two individually tagged aliquots of human PBMCs) and with one CMO ID assigned per sample ID.

Libraries were prepared following the Chromium Next GEM Single Cell 3ʹ Reagent Kits v3.1 (Dual Index) with Feature Barcode technology for Cell Multiplexing User Guide (CG000388) and sequenced on Illumina NovaSeq 6000.

## Download original input data

```bash
curl -O https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/6.0.0/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_Multiplex/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_Multiplex_fastqs.tar
curl -O https://cf.10xgenomics.com/samples/cell-exp/6.0.0/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_Multiplex/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_Multiplex_config.csv
curl -O https://cf.10xgenomics.com/samples/cell-exp/6.0.0/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_Multiplex/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_Multiplex_count_feature_reference.csv
```

## Original `cellranger multi` config

```bash
[gene-expression]
reference,/path/to/references/refdata-gex-GRCh38-2020-A
expect-cells,10000

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_gex,/path/to/fastqs/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_gex,any,Human_PBMC_10K_gex,gene expression,
SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_2_gex,/path/to/fastqs/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_2_gex,any,Human_PBMC_10K_gex,gene expression,
SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_multiplexing_capture,/path/to/fastqs/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_multiplexing_capture,any,Human_PBMC_10K_mult_capture,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
PBMCs_human_1,CMO301,PBMCs_human_1
PBMCs_human_2,CMO302,PBMCs_human_2
```

## Cell Multiplexing Data 

The cell multiplexing libraries require more generous subsampling to yield sufficient reads for successful processing. 
Test data use the first 3,500,000 lines of each of the lane 1 FASTQs.

```bash
[ec2-user@ip-10-0-25-198 cmo]$ zcat ../SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_multiplexing_capture/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_multiplexing_capture_S1_L001_R1_001.fastq.gz | head -n 3500000 | gzip -c > subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_multiplexing_capture_S1_L001_R1_001.fastq.gz
[ec2-user@ip-10-0-25-198 cmo]$ zcat ../SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_multiplexing_capture/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_multiplexing_capture_S1_L001_R2_001.fastq.gz | head -n 3500000 | gzip -c > subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_multiplexing_capture_S1_L001_R2_001.fastq.gz
```
