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

* **mgnifams_input_small.faa**: An amino acid fasta file of metagenomics derived sequences called from the MGnify analysis pipelines. The file contains 50K sequences, a size that allows the pipeline to execute both fast and also create enough clusters/families to sufficiently test all modules of the proteinfamilies pipeline. Called by samplesheets samplesheet.csv and samplesheet_multi_sample_with_gz.csv. It is used in the default, minimal and multi-sample/compressed test configurations.
* **mgnifams_input_small_copy.faa.gz**: A compressed copy of mgnifams_input_small.faa. Called by samplesheet samplesheet_multi_sample_with_gz to simultaneously test for the functionality of both multi-sample and compressed fasta inputs.
* **mgnifams_extra.faa**: An amino acid fasta file of another 50K sequences. Called by samplesheets samplesheet_update.csv and samplesheet_full.csv to test the functionality of the update families mechanism. Sequences that match existing families are processed along those families, which are then updated. Non-hit sequences will go through the basic family generation workflow.
* **existing_hmms.tar.gz**: A compressed archive containing 5 HMM files (.hmm.gz) of previously generated families (from mgnifams_input_small.faa). Called by samplesheets samplesheet_update.csv and samplesheet_full.csv to test the functionality of the update families mechanism.
* **existing_msas.tar.gz**: A compressed archive containing 5 MSA files (.aln) of previously generated families (from mgnifams_input_small.faa). Called by samplesheets samplesheet_update.csv and samplesheet_full.csv to test the functionality of the update families mechanism. The files in the archive are the same in number as those in the HMM archive, and their base file names are matching.
* **mgnifams_extra.faa.gz**: A compressed copy of mgnifams_extra.faa. Called by samplesheets/v3/samplesheet_update.csv to test compressed fasta input in the update families mechanism.

### v3 test data (pipeline v3.0.0 onwards)

Pipeline v3.0.0 changed the samplesheet columns to `id,fasta,existing_hmms,existing_seed_msas,existing_full_msas`. The new samplesheets live in `samplesheets/v3/`; the top-level `samplesheets/` and the archives above are kept unchanged for earlier pipeline releases.

* **samplesheets/v3/samplesheet.csv**, **samplesheet_multi_sample_with_gz.csv**, **samplesheet_full.csv**: Same inputs as their v2 counterparts, with the v3 columns. The `mgnifams_update` sample of samplesheet_full.csv uses the v3 archives below.
* **samplesheets/v3/samplesheet_update.csv**: One sample per update input combination: `update_hmm_only` (HMMs only), `update_hmm_seed_gz` (HMMs + seed MSAs, compressed fasta) and `update_hmm_seed_full` (HMMs + seed MSAs + full MSAs).
* **v3/existing_hmms.tar.gz**: A compressed archive containing 5 HMM files (.hmm.gz). `existing_fam_1` to `existing_fam_4` come from a pipeline v2.6.0 test_full run and all get hits in mgnifams_extra.faa. `existing_fam_zero_hits` is built from 5 sequences of the Pfam PF00087 (snake three-finger toxin) seed and gets no hits in mgnifams_extra.faa, to test the handling of existing families without hits.
* **v3/existing_seed_msas.tar.gz**: A compressed archive containing the 5 seed MSA files (.aln) the HMMs above were built from.
* **v3/existing_full_msas.tar.gz**: A compressed archive containing the 5 full MSA files (.aln) of the same families.

The three v3 archives contain the same families, and their base file names match each other and the `NAME` field of each HMM.

## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).

[^1]: From [stackoverflow](https://stackoverflow.com/a/60846265/11502856)
