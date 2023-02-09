# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the nf-core pipelines

## Introduction

nf-core is a collection of high quality Nextflow pipelines. This repository contains various files for CI and unit testing of nf-core pipelines and infrastructure.

The principle for nf-core test data is as small as possible, as large as necessary. Always ask for guidance on the [nf-core slack](https://nf-co.re/join) before adding new test data.

## Full-size test data

### Assembly test data

The `samplesheet.fullsize_assembly.csv` file links to the following metagenome assembly test data:

| SAMPLE | ASSEMBLY FILE            | CONTIG DEPTHS FILE      |
|--------|--------------------------|-------------------------|
| G0_T0  | fullsize_CAMISIM_SPAdes.G0_T0.top300000.fasta.gz | fullsize_CAMISIM_SPAdes.G0_T0_depths.top300000.tsv |
| G0_T1  | fullsize_CAMISIM_SPAdes.G0_T1.top300000.fasta.gz | fullsize_CAMISIM_SPAdes.G0_T1_depths.top300000.tsv |

The underlying read data was simulated with CAMISIM [(Fritz, A. et al., 2019)](https://doi.org/10.1186/s40168-019-0633-6) based on the genome sources from the "CAMI II challenge toy mouse gut dataset" [(Meyer et al., 2021)](https://doi.org/10.1038/s41596-020-00480-3), containing 791 genomes (the simulated read data is available at https://doi.org/10.5281/zenodo.5155395).
The data was further processed with the nf-core/mag pipeline version 2.1.0. The assemblies computed with `SPADes` as well as the by nf-core/mag generated contig depth files are used here.
Due to the GitHub limit of 100MB, the top 300.000 contigs were extracted.

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
