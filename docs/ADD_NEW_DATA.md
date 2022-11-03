# How to add and use new test dataset

Please fill in the appropriate checklist below (delete whatever is not relevant). These are the most common things requested when adding a new test dataset.

 - [X] Check [here](https://github.com/nf-core/test-datasets/branches/all) that there isn't already a branch containing data that could be used
   - If this is the case, follow the [documentation on how to use an existing test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/USE_EXISTING_DATA.md)
 - [X] Fork the [nf-core/test-datasets repository](https://github.com/nf-core/test-datasets) to your GitHub account
 - [X] Create a new branch on your fork
 - [X] Add your test dataset
   - [ ] If you clone it locally use `git clone <url> --branch <branch> --single-branch`
 - [X] Make a PR on a new branch with a relevant name
 - [ ] Wait for the PR to be merged
 - [ ] Use this newly created branch for your tests

 # Description of dataset

 - genomics/homo_sapiens/genome/genome_config = config file used for cellranger-atac mkref
 - genomics/homo_sapiens/genome/genome_scATAC.gtf = gtf file used for cellranger-atac mkref.  I could not use the gtf file genomics/homo_sapiens/genome/genome.gtf because it seems to contain exons that does belong to transcripts which are not listed in the gtf and thus throwing an error for cellranger-mkref. I therefore filtered the file.
 - genomics/homo_sapiens/genome/genome_motifs.txt = protein motifs used for cellranger-atac mkref

 - genomics/homo_sapiens/illumina/scATAC/test_scATAC_S1_L001_R1_001.fastq.gz = downsampled test data from cellranger-atac page (https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/mkfastq)
 - genomics/homo_sapiens/illumina/scATAC/test_scATAC_S1_L001_R2_001.fastq.gz = downsampled test data from cellranger-atac page (https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/mkfastq)
 - genomics/homo_sapiens/illumina/scATAC/test_scATAC_S1_L001_R3_001.fastq.gz = downsampled test data from cellranger-atac page (https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/mkfastq)
