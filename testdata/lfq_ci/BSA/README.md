### Input spectra (`*.mzML`)
Three runs on a BSA sample. Each in-silico fractionated by splitting the files in RT dimension.
Results in BSA\_[sample]\_[fraction].mzML

### Input database (`*.fasta`)

The 18 proteins supposed to be in each BSA sample. Decoys were already added.
TODO which decoy prefix and which method.

### Experimental design (`testdata/BSA_design.tsv`)

Samples and fractions are listed in the OpenMS specific experimental design format.
