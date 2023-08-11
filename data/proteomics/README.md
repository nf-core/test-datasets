# Test Data

## Table of contents

- [bos_taurus](#bos_taurus)
- [database](#database)
- [maxquant](#maxquant)
- [msspectra](#msspectra)
- [parameter](#parameter)

## bos_taurus

## database

## MaxQuant

For MaxQuant, the files `proteinGroups.txt` and `experimentalDesignTemplate.txt` as well as `SampleKey.docx` were downloaded from the [PXD043349 project in PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD043349).
The contents of `proteinGroups.txt` were left unchanged and the file simply renamed to `MaxQuant_proteinGroups.txt`. The file `experimentalDesignTemplate.txt` was renamed to `MaxQuant_samplesheet.tsv`; also, it was changed like so: A column called Celltype was added; here, for each sample number, its corresponding B cell subset was entered as described in `SampleKey.docx`. Then, as the `proteinGroups.txt` only contained one LFQ intensity column per sample (the Experiment column in `experimentalDesignTemplate.txt`), not per Name (the Name column in `experimentalDesignTemplate.txt`), `MaxQuant_samplesheet.tsv` was reduced to the first row per sample, e.g. AH-5-1 for sample 5, AH-6-1 for sample 6 and so on. This was done as the original purpose of introducing this dataset to nf-core was to use it as test data for the proteus module and to process the LFQ intensity values. Then, the column fakeBatch was added in which the first half of the samples (rounded up, i.e. 8) was classified as b1 and the second half (the last 7 samples) was classified as b2. This is not at all biologically backed (or if it is, only coincidentally) and was simply done to allow for testing batch effect functionalities.
Finally, `MaxQuant_contrasts.csv` was created to be used as a contrast file for the [nf-core/differentialabundance pipeline](https://github.com/nf-core/differentialabundance). Contrasts were decided (for no particular reason) to be set up as comparisons between T1 cells and the other three cell types in the dataset; in the last comparison, the batch effect is also included.

## msspectra

## parameter