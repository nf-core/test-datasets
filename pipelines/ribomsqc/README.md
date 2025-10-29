# ribomsqc test data

## Description
Test data for nf-core/ribomsqc pipeline - ribonucleoside mass spectrometry QC analysis.

## Files
- `testdata/BSA_QC_test.raw`: 2.9MB BSA QC standard mass spectrometry file (.raw format)
- `testdata/samplesheet.csv`: Input samplesheet for ribomsqc pipeline testing
- `testdata/analytes.tsv`: Ribonucleoside analytes definitions for QC analysis

## Usage
Used by ribomsqc pipeline with `-profile test` for automated CI testing:

nextflow run nf-core/ribomsqc -profile test,docker --outdir results

## Data Source
BSA QC standard file from internal dataset of Proteomics Unit CRG for MS analysis validation.
File size optimized for CI testing (2.9MB - fast download/processing).
