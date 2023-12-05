# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

Test data to be used for automated testing with the nf-core pipelines.

This branch contains test data for the [nf-core/crisprseq](https://github.com/nf-core/crisprseq) pipeline.

## Content

### Analysis of CRISPR editing

`testdata-edition` contains the a minimal test dataset needed to run the pipeline for the analysis of CRISPR editing.

- `samplesheet_test.csv` contains the input samplesheet file to run the pipeline.
- The `fastq`files are samples from a simulated run (`chr6`) and from a MiSeq amplicon sequencing run of the target sites `TRAC` and `AAVS1`.

- `samplesheet_test_full.csv` contains the input samplesheet file to run the pipeline with full size data. Data is obtained from ENA project [PRJNA326019](https://www.ebi.ac.uk/ena/browser/view/PRJNA326019) obtained by [VanOverbeek et al. 2016](https://doi.org/10.1016/j.molcel.2016.06.037)

- `samplesheet_test_umis.csv` contains the input samplesheet file to run the pipeline with the option of UMIs. Samples were obtained from ENA project [ERR9897751](https://www.ebi.ac.uk/ena/browser/view/ERR9897751) and [ERR9897765](https://www.ebi.ac.uk/ena/browser/view/ERR9897765) and subsampled to with seqtk and a seed of 100.
  ```
    seqtk sample -s100 ERR9897751_1.fastq 632577 > ERR9897751_1_subsample.fastq
    seqtk sample -s100 ERR9897751_2.fastq 632577 > ERR9897751_2_subsample.fastq
    seqtk sample -s100 ERR9897765_1.fastq 606960 > ERR9897765_1_subsample.fastq
    seqtk sample -s100 ERR9897765_2.fastq 606960 > ERR9897765_2_subsample.fastq

  ```

### Analysis of CRISPR screening

`testdata` contains the a minimal test dataset needed to run the pipeline for the analysis of CRISPR functional screenings.

- `samplesheet_test.csv` contains the input samplesheet file to run the pipeline with single-end data. Data is obtained from ENA project [PRJNA540212](https://www.ebi.ac.uk/ena/browser/view/PRJNA540212) obtained by [Tyner JW et al. 2019](https://doi.org/10.1158/2159-8290.CD-19-0125) and subsampled to with seqtk and a seed of 100.
  ```
    seqtk sample -s100 SRR8983580.fastq.gz 1529806 > SRR8983580.small.fastq.gz
    seqtk sample -s100 SRR8983579.fastq.gz 1529806 > SRR8983579.small.fastq.gz
  ```

- `samplesheet_test_paired.csv` contains the input samplesheet file to run the pipeline with paired-end data. Data is obtained from ENA project [PRJNA924808](https://www.ebi.ac.uk/ena/browser/view/PRJNA924808) obtained by [Yangfan Zhou et al. 2019](https://doi.org/10.1038/s41598-023-38810-6) and subsampled to with seqtk and a seed of 100.
  ```
    seqtk sample -s SEED SRR23101364_1.fastq.gz 0.25 > SRR23101364_1.small.fastq.gz
    seqtk sample -s SEED SRR23101364_2.fastq.gz 0.25 > SRR23101364_2.small.fastq.gz
  ```


- the 'brunello_target_sequence.txt' file is used for the mageck count module and is taken from the mageck [documentation](https://media.addgene.org/cms/filer_public/8b/4c/8b4c89d9-eac1-44b2-bb2f-8fea95672705/broadgpp-brunello-library-contents.txt). Only target transcript,sgRNA Target,Sequence   Target and Gene Symbol are kept as columns.

- 'design_matrix.txt' is the design matrix to run mageck mle and is tailored to the output created by samplesheet_test.csv and mageck count.

- 'count_table.csv' is the count matrix to run MAGeCK MLE. It is taken from [leukemia.new.csv](https://sourceforge.net/projects/mageck/files/example/)

- the folder 'full_test' contains all the necessary files to run the full AWS tests. 
- `samplesheet_full.csv` contains the input samplesheet file to run the pipeline. Data is obtained from ENA project [PRJNA540212](https://www.ebi.ac.uk/ena/browser/view/PRJNA540212) obtained by [Tyner JW et al. 2019](https://doi.org/10.1158/2159-8290.CD-19-0125). 
- the design matrices were manually created.

  ```

