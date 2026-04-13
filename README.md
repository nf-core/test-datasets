# nf-core/denovoproteomics test data

Test data for the [nf-core/denovoproteomics](https://github.com/nf-core/denovoproteomics) pipeline.

## Contents

| File | Description | Size |
|------|-------------|------|
| `samplesheet.csv` | Samplesheet for standard mode (2 samples) | <1 KB |
| `samplesheet_mapping.csv` | Samplesheet for mapping mode (2 samples) | <1 KB |
| `testdata/sample1.mzML` | DDA mzML, 48 spectra (34 MS2), Thermo LTQ FT | 4.9 MB |
| `testdata/sample2.mzML` | DDA mzML, 48 spectra (34 MS2), Thermo LTQ FT | 4.9 MB |

## Cross-branch references

The pipeline test configs also reference data from other nf-core/test-datasets branches:

- **FASTA reference** for mapping mode: `modules` branch, `data/proteomics/database/yeast_UPS_mini.fasta` (4.2 KB, 10 proteins)

## Usage

```bash
# Stub test (CI, no real tools)
nextflow run nf-core/denovoproteomics -profile test -stub --outdir results

# Full test (real tools, small data)
nextflow run nf-core/denovoproteomics -profile test_full --outdir results
```
