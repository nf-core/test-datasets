# Modules test data

## Test profile output data

`nf-core/differentialabundance` was executed with this command:

```bash
nextflow run nf-core/differentialabundance -r 9a84f4bf1e581425c0fcb7edcec772998265eeb1 -latest -profile docker,test --outdir results
```

- `Mus_musculus.anno.feature_metadata.tsv` was extracted from the channel `GTF_TO_TABLE.out.feature_annotation`.

- `Mus_musculus.anno.tsv` was extracted from the channel `VALIDATOR.out.feature_meta`.

- `SRP254919.samplesheet.sample_metadata.tsv` was extracted from the channel `VALIDATOR.out.sample_meta`.

- `all.normalised_counts.tsv` was extracted from the channel `ch_norm`, used as input for the process `CUSTOM_TABULARTOGSEAGCT`.

- `SRP254919.salmon.merged.gene_counts.top1000cov.assay.tsv` was extracted from the channel `ch_matrix_for_differential`, used as input for the process `CUSTOM_MATRIXFILTER`.

- `all.vst.tsv` was extracted from the channel `DESEQ2_DIFFERENTIAL.out.vst_counts`.

- `treatment_mCherry_hND6_.deseq2.results.tsv` was extracted from the channel `DESEQ2_DIFFERENTIAL.out.results`.
  
- `SRP254919.anno.feature_metadata.tsv` was extracted from the channel `GTF_TO_TABLE.out.feature_annotation`.



## Test full profile output data

`nf-core/differentialabundance` was executed with this command:

```bash
nextflow run nf-core/differentialabundance -r  9a84f4bf1e581425c0fcb7edcec772998265eeb1  -profile docker,test_full --outdir results
```

- `study.filtered.tsv` was extracted from the channel `CUSTOM_MATRIXFILTER.out.filtered`.

- `Condition_genotype_WT_KO.deseq2.results_filtered.tsv` was extracted from the channel `CUSTOM_FILTERDIFFERENTIALTABLE.out.filtered`.

- `Condition_genotype_WT_KO.cls` was extracted from the channel `CUSTOM_TABULARTOGSEACLS.out.cls`.

- `Condition_treatment_Control_Treated.gct` was extracted from the channel `CUSTOM_TABULARTOGSEAGCT.out.gct`.

- `Mus_musculus.anno.feature_metadata.chip` was extracted from the channel `TABULAR_TO_GSEA_CHIP.out.chip`.
