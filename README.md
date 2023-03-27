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

### Analysis of CRISPR editing

`testdata` contains the a minimal test dataset needed to run the pipeline for the analysis of CRISPR functional screenings.

- `samplesheet_test.csv` contains the input samplesheet file to run the pipeline.
- The `fastq`files are samples from a simulated run and are downloaded from the [ENA archive](https://www.ebi.ac.uk/ena/browser/view/ERR376998)

- the 'yusa_library.csv' file is used for the mageck count module and is taken from the mageck [documentation](https://sourceforge.net/projects/mageck/files/libraries/)

- 'HT-29_counts.tsv' and 'KY_Library_v1.1.tsv' are both from the [CRISPRcleanR R library](https://github.com/francescojm/CRISPRcleanR/tree/master/data)

- 'design_matrix.txt' is the design matrix to run mageck mle and is tailored to the 'count_table.csv'

- 'count_table.csv' is the count matrix to run mageck mle. It is taken from [leukemia.new.csv](https://sourceforge.net/projects/mageck/files/example/)
