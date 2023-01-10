# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the nf-core pipelines

## Introduction

nf-core is a collection of high quality Nextflow pipelines.
This repository contains various files for CI and unit testing of nf-core pipelines and infrastructure.

The principle for nf-core test data is as small as possible, as large as necessary.
Please see the [guidelines](https://nf-co.re/docs/contributing/test_data_guidelines) for more detailed information. Always ask for guidance on the [nf-core slack](https://nf-co.re/join) before adding new test data.

## Documentation

nf-core/test-datasets comes with documentation in the `docs/` directory:

01. [Add a new  test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/ADD_NEW_DATA.md)
02. [Use an existing test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/USE_EXISTING_DATA.md)

## Downloading nf-core/sarek3 test data

Due the large number of large files in this repository for each pipeline, we highly recommend cloning only the branches you would use.

```bash
git clone https://github.com/nf-core/test-datasets.git --single-branch --branch sarek3
```

## Support

For further information or help, don't hesitate to get in touch on the [Slack `#sarek` channel](https://nfcore.slack.com/channels/sarek) (you can join with [this invite](https://nf-co.re/join/slack)), or contact us: [Gisela Gabernet](mailto:gisela.gabernet@qbic.uni-tuebingen.de?subject=[GitHub]%20nf-core/test-datasets), [Maxime Garcia](mailto:maxime.garcia@seqera.io?subject=[GitHub]%20nf-core/test-datasets), [Friederike Hanssen](mailto:friederike.hanssen@qbic.uni-tuebingen.de?subject=[GitHub]%20nf-core/test-datasets)
