# 10k_pbmc

Peripheral blood mononuclear cells (PBMCs) from a healthy donor were obtained by 10x Genomics from AllCells. 500,000 cells were stained with TotalSeq™-B Human Universal Cocktail, V1.0 (BioLegend, Cat# 399904) following the demonstrated protocol Cell Surface Protein Labeling for Chromium Fixed RNA Profiling for Singleplex Samples with Feature Barcode technology (CG000529, Rev A), with the minor modification that the staining was done in a total of 50μL (BioLegend recommendation for TotalSeq™ Universal Cocktails). After staining, cells were washed using the 2-Wash option and then fixed for 1 hour at room temperature following the demonstrated protocol Fixation of Cells & Nuclei for Chromium Fixed RNA Profiling (CG000478).

Fixed RNA Gene Expression and Cell Surface Protein libraries were generated as described in the Chromium Fixed RNA Profiling for Singleplexed Samples with Feature Barcode technology for Cell Surface Protein User Guide (CG000477) using the Chromium X and sequenced on an Illumina NovaSeq6000 with approximately 37k read pairs per cell.

## Download original input data

```bash
curl -O https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/7.0.0/10k_Human_PBMC_TotalSeqB_singleplex_Multiplex/10k_Human_PBMC_TotalSeqB_singleplex_Multiplex_fastqs.tar
curl -O https://cf.10xgenomics.com/samples/cell-exp/7.0.0/10k_Human_PBMC_TotalSeqB_singleplex_Multiplex/10k_Human_PBMC_TotalSeqB_singleplex_Multiplex_config.csv
curl -O https://cf.10xgenomics.com/samples/cell-exp/7.0.0/10k_Human_PBMC_TotalSeqB_singleplex_Multiplex/10k_Human_PBMC_TotalSeqB_singleplex_Multiplex_count_feature_reference.csv
```

## Original `cellranger multi` config

```bash
[gene-expression]
reference,/path/to/references/refdata-gex-GRCh38-2020-A
expect-cells,10000

[vdj]
reference,/path/to/references/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0

[feature]
reference,/path/to/feature_references/sc5p_v2_hs_PBMC_10k_feature_ref.csv

[libraries]
fastq_id,fastqs,lanes,feature_types,subsample_rate
sc5p_v2_hs_PBMC_10k_5gex,/path/to/fastqs/gex/sc5p_v2_hs_PBMC_10k_5gex_5fb_fastqs/sc5p_v2_hs_PBMC_10k_5gex_fastqs,1|2,gene expression,
sc5p_v2_hs_PBMC_10k_5fb,/path/to/fastqs/gex/sc5p_v2_hs_PBMC_10k_5gex_5fb_fastqs/sc5p_v2_hs_PBMC_10k_5fb_fastqs,1|2,antibody capture,
sc5p_v2_hs_PBMC_10k_b,/path/to/fastqs/vdj/sc5p_v2_hs_PBMC_10k_b_fastqs,1|2,vdj-b,
sc5p_v2_hs_PBMC_10k_t,/path/to/fastqs/vdj/sc5p_v2_hs_PBMC_10k_t_fastqs,1|2,vdj-t,
```
