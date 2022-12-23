# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the nf-core pipelines

> ⚠️ **Do not merge your test data to `master`! Each pipeline has a dedicated branch (and a special one for modules)**

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

## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).

[^1]: From [stackoverflow](https://stackoverflow.com/a/60846265/11502856)

## How data has been acquired
### Initial data
Beware can be long to download especially 'NA12878.final.cram'
```
wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.filtered.shapeit2-duohmm-phased.vcf.gz{,.tbi}
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz hs38DH.chr21.fa.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram
wget https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/chr21/germlineresources/gnomAD.r2.1.1.vcf.gz -O gnomAD.r2.1.1.chr21.vcf.gz
```

### Preparation of the different panel files
```
. get_panel_s.sh
```

### Preparation and downsampling of the individual file validation and test file
```
. get_ind_1x
```

### Compute the genotype likelihood for the individual data 
```
. get_ind_gl.sh
```