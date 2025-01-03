# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the nf-core pipelines

## Introduction

nf-core is a collection of high quality Nextflow pipelines. This repository contains various files for CI and unit testing of nf-core pipelines and infrastructure.

The principle for nf-core test data is as small as possible, as large as necessary. Please see the [guidelines](https://nf-co.re/docs/contributing/test_data_guidelines) for more detailed information. Always ask for guidance on the [nf-core slack](https://nf-co.re/join) before adding new test data.

## Documentation

nf-core/test-datasets comes with documentation in the `docs/` directory:

01. [Add a new  test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/ADD_NEW_DATA.md)
02. [Use an existing test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/USE_EXISTING_DATA.md)

## Downloading test data

Due the large number of large files in this repository for each pipeline, we highly recommend cloning only the branches you would use.

```bash
git clone <url> --single-branch --branch <pipeline/modules/branch_name>
```

To subsequently clone other branches[^1]

```bash
git remote set-branches --add origin [remote-branch]
git fetch
```

## Test data for nf-core/proteinfamilies

This branch contains test data for the [nf-core/proteinfamilies](https://github.com/nf-core/proteinfamilies) pipeline, in the test_data folder.

* **mgnifams_input_small.fa**: An amino acid fasta file of metagenomics derived sequences called from the MGnify analysis pipelines. The file contains 50K sequences, a size that allows the pipeline to execute both fast and also create enough clusters/families to sufficiently test all modules of the proteinfamilies pipeline. Called by samplesheets samplesheet.csv and samplesheet_multi_sample_with_gz.csv. It is used in the default, minimal and multi-sample/compressed test configurations.
* **mgnifams_input_small_copy.fa.gz**: A compressed copy of mgnifams_input_small.fa. Called by samplesheet samplesheet_multi_sample_with_gz to simultaneously test for the functionality of both multi-sample and compressed fasta inputs.
* **mgnifams_extra.fa**: An amino acid fasta file of another 50K sequences. Called by samplesheets samplesheet_update.csv and samplesheet_full.csv to test the functionality of the update families mechanism. Sequences that match existing families are processed along those families, which are then updated. Non-hit sequences will go through the basic family generation workflow.
* **existing_hmms.tar.gz**: A compressed archive containing 5 HMM files (.hmm.gz) of previously generated families (from mgnifams_input_small.fa). Called by samplesheets samplesheet_update.csv and samplesheet_full.csv to test the functionality of the update families mechanism.
* **existing_msas.tar.gz**: A compressed archive containing 5 MSA files (.aln) of previously generated families (from mgnifams_input_small.fa). Called by samplesheets samplesheet_update.csv and samplesheet_full.csv to test the functionality of the update families mechanism. The files in the archive are the same in number as those in the HMM archive, and their base file names are matching.

## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).

[^1]: From [stackoverflow](https://stackoverflow.com/a/60846265/11502856)
