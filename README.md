# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

Test data to be used for automated testing with the nf-core pipelines

> ⚠️ **Do not merge your test data to `master`! Each pipeline has a dedicated branch (and a special one for modules)**

## Introduction

nf-core is a collection of high quality Nextflow pipelines. This repository contains various files for CI and unit testing of nf-core pipelines and infrastructure.

The principle for nf-core test data is as small as possible, as large as necessary. Please see the [guidelines](https://nf-co.re/docs/contributing/test_data_guidelines) for more detailed information. Always ask for guidance on the [nf-core slack](https://nf-co.re/join) before adding new test data.

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

## Datasets for nf-core/proteinannotator

### DIAMOND test data

Created to support the DIAMOND blastp taxonomic classification integration in nf-core/proteinannotator ([PR #50](https://github.com/nf-core/proteinannotator/pull/50)), and to address a CI disk-exhaustion issue flagged during review.

Source: Minimal NCBI taxonomy and RefSeq data subset needed to exercise `DIAMONDPREPARETAXA`, `DIAMOND_MAKEDB`, `DIAMOND_BLASTP`, and `NCBIREFSEQDOWNLOAD` without pulling full production-scale databases. The five accessions below were deliberately chosen from the real `other` category because together they produce genuine DIAMOND alignment.

- `WP_031942563.1` -- tetracycline efflux MFS transporter Tet(B) [Transposon Tn10]
- `WP_430799656.1` -- class D beta-lactamase OXA-1379 [medical waste metagenome]
- `WP_148044478.1` -- phosphoethanolamine--lipid A transferase MCR-5.4 [hospital metagenome]
- `WP_168247882.1` -- extended-spectrum class C beta-lactamase IDC-2 [sediment metagenome]
- `WP_168247881.1` -- extended-spectrum class C beta-lactamase IDC-1 [sediment metagenome]

`test_refseq.fasta`: minimal refseq protein fasta for DIAMOND_BLASTP
`refseq/release/other/other.wp_protein_test.1.protein.faa.gz`: minimal RefSeq 'other' category subset for NCBIREFSEQDOWNLOAD, replacing the full-category download that was exhausting CI disk space
`mini_taxdump.tar.gz`: minimal NCBI taxdump for DIAMONDPREPARETAXA tests
`mini_prot.accession2taxid.gz`: minimal accession2taxid map for DIAMOND_MAKEDB

```
testdata/diamond/mini_taxdump.tar.gz
testdata/diamond/mini_prot.accession2taxid.gz
testdata/diamond/test_refseq.fasta
testdata/diamond/refseq/release/other/other.wp_protein_test.1.protein.faa.gz
```
