# nf-core/denovoproteomics test data

Test data for the [nf-core/denovoproteomics](https://github.com/nf-core/denovoproteomics) pipeline.

## Contents

| File | Description | Size |
|------|-------------|------|
| `samplesheet.csv` | Samplesheet for standard mode (2 samples) | <1 KB |
| `samplesheet_mapping.csv` | Samplesheet for mapping mode (2 samples) | <1 KB |

## Cross-branch references

All spectra files and FASTA references are reused from the `modules` branch to avoid data duplication:

- **Spectra**: `modules` branch, `data/proteomics/msspectra/OVEMB150205_12.raw` (22.5 MB) and `OVEMB150205_14.raw` (26.5 MB)
- **FASTA reference** for mapping mode: `modules` branch, `data/proteomics/database/yeast_UPS_mini.fasta` (4.2 KB, 10 proteins)

## Usage

```bash
# Stub test (CI, no real tools)
nextflow run nf-core/denovoproteomics -profile test -stub --outdir results

# Full test (real tools, small data)
nextflow run nf-core/denovoproteomics -profile test_full --outdir results
```
