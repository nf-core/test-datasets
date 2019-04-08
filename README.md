# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the nf-core pipelines

## Testdataset for denovohybrid pipeline

The testdata folder contains the following files that are used with the [denovohybrid](https://github.com/nf-core/denovohybrid) pipeline.

Subset of illumina read pairs from Bacterial sample

- `test_read_illumina_R1.fastq.gz` 
- `test_read_illumina_R2.fastq.gz` 

Matching subset of Nanopore reads from same sample

- `test_read_nanopore.fastq.tz`

Input files with read paths in .tsv format. These files are provided to the pipeline with the `--input` parameter. 

- `test_files_lr.tsv`
- `test_files.tsv`

One file is used for the nanopore only test, the other for hybrid assembly

## Support

For further information or help, don't hesitate to get in touch on our [Gitter channel](https://gitter.im/nf-core/Lobby)
