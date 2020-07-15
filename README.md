# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the nf-core pipelines

## Introduction

This is the gwas example-data branch, part of the nf-core collection of high quality Nextflow pipelines.

## Git clone the gwas pipeline test data

If you want to get a local copy of the test data, you can either git clone the whole test data material, including all test data for all nf-core pipelnies, or if you want to save storage space you can clone the example data for one specific pipeline. 

The data in this example-data branch is the same as the gwas pipeline uses for testing. It is accessed simply by cloning the branch either directly from nf-core if you just want to access the data, or if you want to update the data and make pull-request, it is suggested that you first fork the repository and then clone from your personal fork.

```
# If you are a normal user that wants to get a local copy of the test data
git clone -b gwas --single-branch git@github.com:nf-core/test-datasets.git

# If you are a developer and want to update the test data, fork first and then
#  use this command, substituting with your github username
git clone -b gwas --single-branch git@github.com:USERNAME/test-datasets.git

```

## Documentation

nf-core/test-datasets comes with documentation in the `docs/` directory and scripts to generate the example data in the `scripts/` directory.

## Example data organisation
nf-core/test-datasets generated test data is located in the `data/` directory.

```
.
├── data_phenotypes_and_covariates
│   ├── example1.covar
│   └── example1.pheno
├── data_shrink_chunk_4500
│   ├── chr10.vcf.bgz
│   ├── chr10.vcf.bgz.tbi
│   ├── chr11.vcf.bgz
│   ├── chr11.vcf.bgz.tbi
│
└── data_shrink_combined_4500
    ├── chr1_to_22_and_X.vcf.bgz
    └── chr1_to_22_and_X.vcf.bgz.tbi
```

## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).

