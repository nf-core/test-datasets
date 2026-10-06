# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

Test data to be used for automated testing with the nf-core pipelines

> ⚠️ **Do not merge your test data to `master`! Each pipeline has a dedicated branch (and a special one for modules)**

## Introduction

nf-core is a collection of high quality Nextflow pipelines. This repository contains various files for CI and unit testing of nf-core pipelines and infrastructure.

The principle for nf-core test data is as small as possible, as large as necessary. Please see the [guidelines](https://nf-co.re/docs/contributing/test_data_guidelines) for more detailed information. Always ask for guidance on the [nf-core slack](https://nf-co.re/join) before adding new test data.

## nf-core/deepmutscan test data

| File                                                          | Description                                                                                                                           |
| ------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| `testdata/GID1A.fasta`                                        | Reference amplicon containing the GID1A ORF (UniProt Q9MAA7); the mutagenised reading frame is `352-1383`.                            |
| `testdata/reads/GID1A_{input,output}{1,2}_50k_{1,2}.fastq.gz` | 50,000 read-pair subsamples of the GluePCA GID1A-GAI shotgun DMS libraries below, used by `-profile test`.                            |
| `samplesheet/GID1A_test.csv`                                  | Samplesheet for `-profile test` (2 input + 2 output replicates, the subsampled reads).                                                |
| `samplesheet/GID1A_full.csv`                                  | Samplesheet for `-profile test_full`: all 6 libraries (3 input, 3 output at 2500 uM GA3), ~480M read pairs, linked directly from ENA. |

Source: ENA project [PRJEB110196](https://www.ebi.ac.uk/ena/browser/view/PRJEB110196) (runs ERR16945046-ERR16945051), publicly available under the INSDC data policy. The full samplesheet points at the submitted FASTQ files, which ENA archives for all six runs.

## Documentation

nf-core/test-datasets comes with documentation in the `docs/` directory:

1.  [Add a new test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/ADD_NEW_DATA.md)
2.  [Use an existing test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/USE_EXISTING_DATA.md)

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
