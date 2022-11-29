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

| SRA ID | Sex        | Population | Sample title                                                                  |
|---------------|------------|------------|-------------------------------------------------------------------------------|
| ERR188383     | male       | GBR        | ChrX subsampled total RNA-seq male GBR (British from England)                 |
| ERR188428     | female     | GBR        | ChrX subsampled total RNA-seq female GBR (British from England)               |
| ERR188454     | male       | YRI        | ChrX subsampled total RNA-seq male YRI (Yoruba from Ibadan, Nigeria)          |
| ERR204916     | female     | YRI        | ChrX subsampled total RNA-seq female YRI (Yoruba from Ibadan, Nigeria) |

### Sampling procedure

ChrX fastq files we further sub-sampled to improve speed of testing.

1. The example command below was used to sub-sample the raw paired-end FastQ files to 50,000 reads (see [seqtk](https://github.com/lh3/seqtk)):

  ```console
  seqtk sample -s100 ERR188383_unclass_chrX_1.fastq.gz 50000 | gzip > ERR188383_chrX_1.fastq.gz
  seqtk sample -s100 ERR188383_unclass_chrX_2.fastq.gz 50000 | gzip > ERR188383_chrX_2.fastq.gz
  ```
 2. Ensembl GRCh37 annotation was downloaded from iGenomes (https://ewels.github.io/AWS-iGenomes/) and subsetted for ChrX with the following command:
 
  ```console
  awk -F"\t" '$1=="X"' genes.gtf > genes_chrX.gtf
  ```
 
## Full test dataset origin

*H. sapiens* total RNA paired-end RNA-seq dataset was obtained from:

> Shen S, Park JW, Lu ZX, Lin L, Henry MD, Wu YN, Zhou Q, Xing Y. rMATS: robust and flexible detection of differential alternative splicing from replicate RNA-Seq data. Proc Natl Acad Sci U S A. 2014 Dec 23;111(51):E5593-601. doi: 10.1073/pnas.1419161111. [Pubmed](https://pubmed.ncbi.nlm.nih.gov/25480548/) [SRA](https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP014759)

### Sample information full test dataset

| SRA ID    | Title           | Platform            | Sample title                                                        |
|-----------|-----------------|---------------------|---------------------------------------------------------------------|
| SRR536350 | PC3E sample     | Illumina HiSeq 2000 | RNA-seq of prostate cancer cell line E-cadherin positive PC-3 cells |
| SRR536348 | PC3E sample     | Illumina HiSeq 2000 | RNA-seq of prostate cancer cell line E-cadherin positive PC-3 cells |
| SRR536342 | PC3E sample     | Illumina HiSeq 2000 | RNA-seq of prostate cancer cell line E-cadherin positive PC-3 cells |
| SRR536346 | GS689.Li sample | Illumina HiSeq 2000 | RNA-seq of prostate cancer cell line isolated from a secondary metastatic liver tumor after intravenous injection of PC-3 cells into mouse |
| SRR536352 | GS689.Li sample | Illumina HiSeq 2000 | RNA-seq of prostate cancer cell line isolated from a secondary metastatic liver tumor after intravenous injection of PC-3 cells into mouse |
| SRR536344 | GS689.Li sample | Illumina HiSeq 2000 | RNA-seq of prostate cancer cell line isolated from a secondary metastatic liver tumor after intravenous injection of PC-3 cells into mouse |
