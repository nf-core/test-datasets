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

## nf-core/variantbenchmarking test data

samplesheet/

samplesheet_small_germline_hg38: Sample sheet for nf-core/variantbenchmarking test profiles for hg38 small (snvs and indels) variants from germline sample HG002
samplesheet_sv_germline_hg38: Sample sheet for nf-core/variantbenchmarking test profiles for hg38 structural variants from germline sample HG002
samplesheet_indel_somatic_hg38: Sample sheet for nf-core/variantbenchmarking test profiles for hg38 indel variants from somatic sample SEQC2
samplesheet_snv_somatic_hg38: Sample sheet for nf-core/variantbenchmarking test profiles for hg38 snv variants from somatic sample SEQC2
samplesheet_sv_somatic_hg38: Sample sheet for nf-core/variantbenchmarking test profiles for hg38 stuctural variants from somatic sample SEQC2
samplesheet_sv_somatic_hg37_liftover: Sample sheet for nf-core/variantbenchmarking test profiles for hg37 stuctural variants from HG002 sample for lifting over variants


## test data

### germline 

HG002 GiAB sample is used for germline benchmarking

-- test case

*hg37* 

- delly, lumpy and manta are downloaded from https://zenodo.org/records/10428664 

- svaba calls are downloaded from https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/Broad_svaba_05052017/ 

*hg38*

- Ashkenaizm result is from https://github.com/CenterForMedicalGeneticsGhent using HG002

- dragen sample is from https://zenodo.org/records/10428664 

- manta sample is from https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/BU_GRCh38_SVs_06252018/manta.HG002.vcf.gz

- lumpy sample is from https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/BU_GRCh38_SVs_06252018/ajtrio.lumpy.svtyper.HG002.md.sorted.recal.vcf.gz


### somatic

SEQC2 is used for germline benchmarking

*hg38*

- HCC1395T_vs_HC1395N analysis is the results of nf-core/sarek (v3.4.2) with hg38 for SEQC2 tests.

- SEQC2 CNV test files are downloaded from https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/analysis/CNVs/WES/

## truth data

*hg37* 

- data downloaded from https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/ 

*hg38* 

- HG002 data downloaded from https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_HG002_medical_genes_SV_benchmark_v0.01/ 

- SEQC2 data is downloaded from https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/

- SEQC2 CNV benchmarks are downloaded from https://zenodo.org/records/14619054

!Note that all the original files are downsized to chromosome 21 for test purposes. 
