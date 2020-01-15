# test-datasets `proteomicslfq`
Test data to be used for automated testing with the nf-core pipeline proteomicslfq

## Content of this repository

`testdata`, contains all inputs needed for a basic test of proteomicslfq.

### Input spectra (`testdata/*.mzML`)
Three runs on a BSA sample. Each in-silico fractionated by splitting the files in RT dimension.
Results in BSA\_[sample]\_[fraction].mzML

### Input database (`testdata/*.fasta`)

The 18 proteins supposed to be in each BSA sample. Decoys were already added.
TODO which decoy prefix and which method.

### Experimental design (`testdata/BSA_design.tsv`)

Samples and fractions are listed in the OpenMS specific experimental design format.

