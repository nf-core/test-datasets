# File Generation Documentation

This document explains how specific files in the folder were generated during the execution of the **nf-core/nanostring** pipeline.

## Pipeline Details

- **Test Profile**: `test`
- **Revision**: `8efb6b1412805cbd056aee16b6cf504c6a1b716a`

The pipeline was executed with the following command:

```
nextflow run nf-core/nanostring -r 8efb6b1412805cbd056aee16b6cf504c6a1b716a -profile docker,test --outdir results
```

## Generated Files

1. **`normalized_counts.tsv`**
   - Extracted from the outputs of the process `NACHO_NORMALIZE`.

2. **`counts_Norm_GEX_ENDO.tsv`**
   - Extracted from the outputs of the process `CREATE_ANNOTATED_TABLES`.
