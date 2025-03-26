# Test Data

## Table of contents

- [database](#database)
- [maxquant](#maxquant)
- [msspectra](#msspectra)
- [parameter](#parameter)
- [pdb](#pdb)

## database
'UP000005640_9606.fasta' is the reviewed human proteome of the [SWISS-PROT](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102476/) and was downloaded from UniProt.

## maxquant

For MaxQuant, the files `proteinGroups.txt` and `experimentalDesignTemplate.txt` as well as `SampleKey.docx` were downloaded from the [PXD043349 project in PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD043349).
The contents of `proteinGroups.txt` were left unchanged and the file simply renamed to `MaxQuant_proteinGroups.txt`. The file `experimentalDesignTemplate.txt` was renamed to `MaxQuant_samplesheet.tsv`; also, it was changed like so: A column called Celltype was added; here, for each sample number, its corresponding B cell subset was entered as described in `SampleKey.docx`. Then, as the `proteinGroups.txt` only contained one LFQ intensity column per sample (the Experiment column in `experimentalDesignTemplate.txt`), not per Name (the Name column in `experimentalDesignTemplate.txt`), `MaxQuant_samplesheet.tsv` was reduced to the first row per sample, e.g. AH-5-1 for sample 5, AH-6-1 for sample 6 and so on. This was done as the original purpose of introducing this dataset to nf-core was to use it as test data for the proteus module and to process the LFQ intensity values. Then, the column fakeBatch was added in which the first half of the samples (rounded up, i.e. 8) was classified as b1 and the second half (the last 7 samples) was classified as b2. This is not at all biologically backed (or if it is, only coincidentally) and was simply done to allow for testing batch effect functionalities.
Finally, `MaxQuant_contrasts.csv` was created to be used as a contrast file for the [nf-core/differentialabundance pipeline](https://github.com/nf-core/differentialabundance). Contrasts were decided (for no particular reason) to be set up as comparisons between T1 cells and the other three cell types in the dataset; in the last comparison, the batch effect is also included. Additionally, in the last row, a contrast was created out of the fakeBatch column; this was done simply to allow for module/pipeline tests with multiple conditions.

To test some downstream modules, the proteus R package was used to read the proteinGroups and create an abundance matrix. For this, on a MacBook Pro with a 2,4 GHz Quad-Core Intel Core i5 processor, in RStudio version 2023.06.2+561 and R version 4.2.3 (2023-03-15), the following script was run using Proteus version 0.2.16 in order to produce `proteus.raw_MaxQuant_proteingroups_tab.tsv`:

'''
library(proteus)
sample.sheet <- read.table("MaxQuant_samplesheet.tsv", header=T, sep="\t") # https://github.com/nf-core/test-datasets/raw/modules/data/proteomics/maxquant/MaxQuant_samplesheet.tsv
sample.sheet$sample <- sample.sheet[["Experiment"]]
sample.sheet$condition <- sample.sheet[["Celltype"]]
measure.cols <- setNames(paste0("LFQ intensity ", sample.sheet[["Experiment"]]), sample.sheet[["Experiment"]])
proteinGroups <- readProteinGroups(file="MaxQuant_proteinGroups.txt", meta=sample.sheet, measure.cols=measure.cols, data.cols=proteus::proteinColumns) # https://github.com/nf-core/test-datasets/raw/modules/data/proteomics/maxquant/MaxQuant_proteinGroups.txt
proteinGroups$tab <- round(log2(proteinGroups$tab), digits=3)
out_df <- data.frame(proteinGroups$tab, check.names = FALSE)
out_df[["Majority protein IDs"]] <- rownames(proteinGroups$tab)
out_df <- out_df[c("Majority protein IDs", colnames(out_df)[colnames(out_df) != "Majority protein IDs"])]
write.table(out_df, file = 'proteus.raw_MaxQuant_proteingroups_tab.tsv', row.names = FALSE, sep = '\t', quote = FALSE) # https://github.com/nf-core/test-datasets/raw/modules/data/proteomics/maxquant/proteus.raw_MaxQuant_proteingroups_tab.tsv
'''

## msspectra
'PXD012083_e005640_II.raw' is a Thermo Fisher RAW file downloaded from [PRIDE](https://www.ebi.ac.uk/pride/) using the project ID PXD012083.
'peakpicker_tutorial_1.mzML' is a mass spectrum file in the open mzML format. The file got retrieved from the [OpenMS](https://github.com/OpenMS/OpenMS) test data on GitHub

## parameter

## pdb

The pdb folder contains protein structure files in .PDB format. Files 1tim.pdb and 8tim.pdb are part of the example datasets used in the foldseek tool (https://github.com/steineggerlab/foldseek).
They describe chicken muscle proteins (engineered and breast respectively) and their structures were determined through X-ray diffraction.
