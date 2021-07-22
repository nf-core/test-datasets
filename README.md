# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the nf-core pipelines

## Introduction

nf-core is a collection of high quality Nextflow pipelines.

## Documentation
nf-core/test-datasets comes with documentation in the `docs/` directory:

01. [Add a new  test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/ADD_NEW_DATA.md)
02. [Use an existing test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/USE_EXISTING_DATA.md)

## Single-cell WGS test dataset

The test dataset is taken from [DOI:10.1002/smll.202001172](https://doi.org/10.1002/smll.202001172). The example command below was used to subsample the raw paired-end FastQ files to 30,000 reads:

```
seqtk sample -s100 A10_combined_R1.fastq.gz 30000 | gzip > A10_test_R1.fastq.gz
seqtk sample -s100 A10_combined_R2.fastq.gz 30000 | gzip > A10_test_R2.fastq.gz
```

## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).
