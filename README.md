# test-datasets: `rnasplice`

This branch contains test data to be used for automated testing with the [nf-core/rnasplice](https://github.com/nf-core/rnasplice) pipeline.

## Content of this repository

`samplesheet/samplesheet.csv`: Experiment design file for minimal test dataset

`samplesheet/samplesheet_full.csv`: Experiment design file for full test dataset

`reference/`: Sub-sampled genome reference files (ChrX)

`testdata/*.fastq.gz`: Subsampled fastq files (ChrX)

## Minimal test dataset origin

ChrX subsampled *H.sapiens* total RNA paired-end RNA-seq data was obtained from:

> Pertea, M., Kim, D., Pertea, G. et al. Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. Nat Protoc 11, 1650–1667 (2016). doi: 10.1038/nprot.2016.095. [Pubmed](https://pubmed.ncbi.nlm.nih.gov/27560171/)

Original data source: 

> Lappalainen T, Sammeth M, Friedländer MR, 't Hoen PA, Monlong J, Rivas MA, Gonzàlez-Porta M, Kurbatova N, Griebel T, Ferreira PG, Barann M, Wieland T, Greger L, van Iterson M, Almlöf J, Ribeca P, Pulyakhina I, Esser D, Giger T, Tikhonov A, Sultan M, Bertier G, MacArthur DG, Lek M, Lizano E, Buermans HP, Padioleau I, Schwarzmayr T, Karlberg O, Ongen H, Kilpinen H, Beltran S, Gut M, Kahlem K, Amstislavskiy V, Stegle O, Pirinen M, Montgomery SB, Donnelly P, McCarthy MI, Flicek P, Strom TM; Geuvadis Consortium, Lehrach H, Schreiber S, Sudbrak R, Carracedo A, Antonarakis SE, Häsler R, Syvänen AC, van Ommen GJ, Brazma A, Meitinger T, Rosenstiel P, Guigó R, Gut IG, Estivill X, Dermitzakis ET. Transcriptome and genome sequencing uncovers functional variation in humans. Nature. 2013 Sep 26;501(7468):506-11. doi: 10.1038/nature12531. [Pubmed](https://pubmed.ncbi.nlm.nih.gov/24037378/) [SRA](https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=study&acc=ERP001942)

Of all samples only the following were chosen for testing: ERR188383, ERR188428, ERR188454 and ERR204916.

### Sampling information

| run_accession | sex        | population | sample_title                                                                  |
|---------------|------------|------------|-------------------------------------------------------------------------------|
| ERR188383     | male       | GBR        | ChrX subsampled total RNA-seq male GBR (British from England)                 |
| ERR188428     | female     | GBR        | ChrX subsampled total RNA-seq female GBR (British from England)               |
| ERR188454     | male       | YRI        | ChrX subsampled total RNA-seq male YRI (Yoruba from Ibadan, Nigeria)          |
| ERR204916     | female     | YRI        | ChrX subsampled total RNA-seq from a female YRI (Yoruba from Ibadan, Nigeria) |

### Sampling procedure

Although chrX had already been sampled from these fastq files we also further sub-sampled them for improved speed of testing.

1. The example command below was used to sub-sample the raw paired-end FastQ files to 50,000 reads (see [seqtk](https://github.com/lh3/seqtk)):

  ```console
  seqtk sample -s100 ERR188383_unclass_chrX_1.fastq.gz 50000 | gzip > ERR188383_chrX_1.fastq.gz
  seqtk sample -s100 ERR188383_unclass_chrX_2.fastq.gz 50000 | gzip > ERR188383_chrX_2.fastq.gz
  ```
