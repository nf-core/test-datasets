# nf-core/denovoproteomics test data

Test data for the [nf-core/denovoproteomics](https://github.com/nf-core/denovoproteomics) pipeline.

## Contents

| File | Description | Size |
|------|-------------|------|
| `samplesheet.csv` | Samplesheet for standard mode (2 samples) | <1 KB |
| `samplesheet_mapping.csv` | Samplesheet for mapping mode (2 samples) | <1 KB |
| `vendor/bruker_timstof_dia.d/` | Bruker timsTOF diaPASEF acquisition (TDF) | 11 MB |
| `vendor/sciex_qtrap.wiff` + `.wiff.scan` | Sciex QTRAP acquisition | 3.3 MB |
| `winnow/winnow_psms.csv` | 500 winnow-scored PSMs, the input to protein assembly | 57 KB |

## Vendor spectra

These exercise the conversion modules, which are the only part of the pipeline
that touches vendor binary formats. Conversion is independent of acquisition
mode, so a DIA acquisition tests a converter just as well as a DDA one.

### `vendor/bruker_timstof_dia.d`

Drives `TDF2MZML`. 500 frames, TDF schema 3; converts to an mzML carrying
30 MS1 and 940 MS2 spectra.

**Source:** [tacular-omics/tdfextractor](https://github.com/tacular-omics/tdfextractor),
`tests/data/example_dia.d`, MIT licence, Copyright (c) 2023 Patrick Garrett.

A `.d` only converts if it holds TDF data with a populated `GlobalMetadata`
table. The MannLabs/timsrust fixtures are far smaller but are simulated and lack
that table, so tdf2mzml rejects them; the other small `.d` archives in
circulation are BAF (QTOF), which is a different format entirely.

A DDA alternative, `200ngHeLaPASEF_1min.d`, is available from the same source at
60 MB. It converts cleanly to 65 MS1 and 2519 MS2 spectra, but at 50 MB of mzML
it is too slow for CPU-only de novo prediction to finish in a test.

### `vendor/sciex_qtrap.wiff` + `vendor/sciex_qtrap.wiff.scan`

Drives `MSCONVERT`. Converts to a 9.9 MB mzML carrying 2127 MS1 and 108 MS2
spectra. A `.wiff` is only readable alongside its `.wiff.scan` companion, so
both files are needed and the module stages them together.

**Source:** [ProteoWizard/pwiz](https://github.com/ProteoWizard/pwiz),
`pwiz_tools/BiblioSpec/tests/inputs/201208-378803.wiff`, Apache-2.0.

Note that msconvert names its output after the sample recorded inside the
bundle (`sciex_qtrap-ABRR-AUG-1.mzML` here), not after the input file, and a
`.wiff` holding several samples yields one mzML per sample.

The corresponding test is excluded from CI: msconvert runs the closed-source
Sciex reader under wine in a 6.7 GB image carrying a vendor licence agreement
that has to be accepted by hand.

## Assembly input

### `winnow/winnow_psms.csv`

500 winnow-scored PSMs with the columns protein assembly consumes:
`spectrum_id`, `prediction`, `calibrated_confidence`, `psm_fdr`, `psm_q_value`,
`psm_pep`. Lets the assembly subworkflow be tested without first running
prediction and rescoring.

## Cross-branch references

Spectra and FASTA references not listed above are reused from the `modules`
branch to avoid data duplication:

- **Spectra**: `data/proteomics/msspectra/OVEMB150205_12.raw` (22.5 MB) and
  `OVEMB150205_14.raw` (26.5 MB)
- **FASTA reference** for mapping mode: `data/proteomics/database/yeast_UPS_mini.fasta`
  (4.2 KB, 10 proteins)

## Usage

```bash
# Stub test (CI, no real tools)
nextflow run nf-core/denovoproteomics -profile test -stub --outdir results

# Full test (real tools, small data)
nextflow run nf-core/denovoproteomics -profile test_full --outdir results
```
