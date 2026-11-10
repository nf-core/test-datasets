# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

Test data to be used for automated testing with the nf-core pipelines

> ⚠️ **Do not merge your test data to `master`! Each pipeline has a dedicated branch (and a special one for modules)**

## Introduction

This branch contains test data for the [nf-core/genephylo](https://github.com/nf-core/genephylo) pipeline.

## Documentation

nf-core/test-datasets comes with documentation in the `docs/` directory:

1.  [Add a new test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/ADD_NEW_DATA.md)
2.  [Use an existing test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/USE_EXISTING_DATA.md)

## Downloading test data

Due the large number of large files in this repository for each pipeline, we highly recommend cloning only the branch you would use.

```bash
git clone https://github.com/nf-core/test-datasets.git --single-branch --branch genephylo
```

To subsequently clone other branches[^1]

```bash
git remote set-branches --add origin [remote-branch]
git fetch
```

## Content of this repository

`testdata/blastdb/mammals_db.fasta`: 60 full-length mRNAs (5'UTR + CDS + 3'UTR, 1.7 to 2.7 kb) from 15 mammal species, used to build a small nucleotide BLAST database.

`testdata/blastdb/mammals_taxid_map.txt`: two-column, space-separated mapping `accession.version taxid` for every sequence in `mammals_db.fasta`, for use with `makeblastdb -taxid_map`.

`testdata/blastquery/TSTG1_nt.fasta`, `testdata/blastquery/TSTG2_nt.fasta`: CDS (nucleotide) of human TSTG1 and TSTG2 respectively. Queries for `blastn`.

`testdata/blastquery/TSTG1_aa.fasta`, `testdata/blastquery/TSTG2_aa.fasta`: translated CDS (protein) of human TSTG1 and TSTG2 respectively. Queries for `tblastn`.

`testdata/speciestree/mammals_sptree.nwk`: rooted species tree (Newick, topology only, no branch lengths) of the 15 species in `mammals_db.fasta`, with full species names as leaves (spaces replaced by underscores, e.g. `Homo_sapiens`). For reconciliation with `treerecs`.

`testdata/speciestree/mammals_sptree_codes.nwk`: same tree with species codes as leaves (first 3 letters of genus + first 3 letters of species, e.g. `HOMSAP`; for three-part names, 3 + 2 + 1 letters, e.g. `CANLUF` for _Canis lupus familiaris_).

### Species

| Species               | NCBI taxid | Sequences |
| --------------------- | ---------- | --------- |
| Homo sapiens          | 9606       | 5         |
| Pan troglodytes       | 9598       | 4         |
| Gorilla gorilla       | 9595       | 4         |
| Macaca mulatta        | 9544       | 4         |
| Callithrix jacchus    | 9483       | 3         |
| Mus musculus          | 10090      | 5         |
| Rattus norvegicus     | 10116      | 3         |
| Oryctolagus cuniculus | 9986       | 4         |
| Bos taurus            | 9913       | 5         |
| Ovis aries            | 9940       | 3         |
| Sus scrofa            | 9823       | 4         |
| Equus caballus        | 9796       | 4         |
| Canis lupus           | 9615       | 5         |
| Felis catus           | 9685       | 3         |
| Loxodonta africana    | 9785       | 4         |

### Genes

| Gene  | Description                              | Species |
| ----- | ---------------------------------------- | ------- |
| TSTG1 | test gene family member 1                | 15      |
| TSTG2 | test gene family member 2 (paralog)      | 13      |
| TSTG3 | test gene family member 3 (distant)      | 12      |
| TSTD1 | test decoy gene 1 (unrelated to TSTG1-3) | 10      |
| TSTD2 | test decoy gene 2 (unrelated to TSTG1-3) | 10      |

## Test dataset origin

The sequences are **simulated** and are not real biological sequences.

Accessions follow RefSeq style (`NM_` for human, mouse, rat and cow; `XM_` "PREDICTED" records for the other species) but are made up and may collide with real RefSeq records. FASTA headers follow the RefSeq format, e.g.:

```text
>XM_062602053.1 PREDICTED: Callithrix jacchus test gene family member 3 (TSTG3), mRNA
```

### Expected results

Verified with BLAST+ 2.12, using human TSTG1 or TSTG2 as query and `-task blastn`:

- `-perc_identity 80`: hits TSTG1 and TSTG2 in 13 species, no TSTG3 (paralog hits at 86 to 89% identity).
- `-perc_identity 70`: additionally hits TSTG3 in 12 species (hits at 73.5 to 76% identity).
- The decoys (TSTD1, TSTD2) never hit.

## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).

[^1]: From [stackoverflow](https://stackoverflow.com/a/60846265/11502856)
