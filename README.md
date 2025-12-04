# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the nf-core pipelines

## Test data and references variantprioritization
This branch contains test data and references for the [nf-core/variantprioritization](https://github.com/nf-core/variantprioritization) pipeline.

## Content of this repository

`test_data/*`: This folder contains minimal vcf files for chr22 with their index and a cna file. 

`samplesheet/default.csv`: Experiment design file for minimal test dataset.

`reference`:
    - `vep_cache_113_GRCh38_chr22.tar.gz`: VEP Cache downsampled to chr22 and with only 1% of all entries in `all_vars.gz` kept for CI testing.
    - `subsample_all_vars.sh`: Script to reduce entries in all_vars.gz file (for VEP cache)
    - `index_subsample.sh`: Index and and zip the subsample all_vars file (for VEP cache)
    - `pcgr_ref_grch38_chr22.tar.gz`: PCGR database dump from `20250314` downsampled to chr22 and with only 0.05% of all entries in all tables kept for CI testing
    - `downsample_pcgr.sh`: Script to reduce tables to 0.05% (for PCGR database)

### Downsampling VEP Cache
```bash
curl -O https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz
tar xzf homo_sapiens_vep_113_GRCh38.tar.gz

mkdir vep_cache_113_GRCh38_chr22
mkdir vep_cache_113_GRCh38_chr22/homo_sapiens
mkdir vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/

cp -r homo_sapiens/113_GRCh38/MT vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/
cp -r homo_sapiens/113_GRCh38/22 vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/
cp homo_sapiens/113_GRCh38/chr_synonyms.txt vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/
cp homo_sapiens/113_GRCh38/info.txt vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/

bash subsample_all_vars.sh
# mamba activate bcf (next script needs bgzip)
bash index_subsample.sh

rm -rf vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/22/subsampled_vars
rm -rf vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/22/all_vars.gz
rm -rf vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/22/all_vars.gz.csi

mv vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/22/subsampled_vars.sorted.bgz vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/22/all_vars.gz
mv vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/22/subsampled_vars.sorted.bgz.csi vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/22/all_vars.gz.csi

tar cvf - vep_cache_113_GRCh38_chr22 | gzip -v > vep_cache_113_GRCh38_chr22.tar.gz
```

### Downsample PCGR database
```bash
BUNDLE_VERSION="20250314"
GENOME="grch38"
BUNDLE="pcgr_ref_data.${BUNDLE_VERSION}.${GENOME}.tgz"
wget https://insilico.hpc.uio.no/pcgr/${BUNDLE}
gzip -dc ${BUNDLE} | tar xvf -

mkdir ${BUNDLE_VERSION}
mv data/ ${BUNDLE_VERSION}

mkdir chr22/
mkdir chr22/20250314/
mkdir chr22/20250314/data/
mkdir chr22/20250314/data/grch38

bash downsample_pcgr.sh

cd chr22
tar cvf - 20250314 | gzip -v > pcgr_ref_grch38_chr22.tar.gz
```

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
