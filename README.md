# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the nf-core pipelines

## Introduction

nf-core is a collection of high quality Nextflow pipelines. This repository contains various files for CI and unit testing of nf-core pipelines and infrastructure.

The principle for nf-core test data is as small as possible, as large as necessary. Always ask for guidance on the [nf-core slack](https://nf-co.re/join) before adding new test data.

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

## Datasets for nf-core/phyloplace

### `testdata/hmmrank/` - HMMER search tables for the `hmmer/hmmrank` module tests

Four gzipped HMMER 3.4 output tables, two per profile, for the `NrdJ` and `NrdJm1` ribonucleotide reductase profiles:

- `NrdJ.tbl.gz`, `NrdJm1.tbl.gz`: per-sequence hit tables (`--tblout`)
- `NrdJ.domtbl.gz`, `NrdJm1.domtbl.gz`: per-domain hit tables (`--domtblout`)

They are excerpts of a real `hmmsearch` of both profiles against protein sequences predicted from the GTDB r226 representative genomes, run through nf-core/phyloplace in "search and place" mode.
Five target sequences were selected and every row for them kept verbatim, giving nine `(target, profile)` pairs, four of the five targets being hit by both profiles.
The trailing comment block that `hmmsearch` writes at the end of each table was removed, because it records the run's working directory and a timestamp and so is neither reproducible nor meaningful outside the original run.
Nothing else was altered.

The five targets were picked to cover the cases that the per-domain coordinate arithmetic has to get right: a pair matching in a single domain, a pair with several overlapping domains, a pair with several disjoint domains, and a pair whose domains overlap in profile coordinates while lying far apart on the target sequence.
That last case is the point of the selection: `AACBMDAH_04658` against `NrdJm1` spans positions 14 to 803 on the sequence but covers only 400 of them, which is what a gene split over more than one open reading frame looks like in these tables.

Consumed by the `hmmer/hmmrank` module tests in nf-core/modules.

### `testdata/PF14720_seed_embedded_taxonomy.alnfaa` - embedded-taxonomy variant of `PF14720_seed.alnfaa`

Same 138 sequences and alignment as `PF14720_seed.alnfaa` (the reference alignment several other phyloplace test profiles already use), with each header's existing taxonomy string from `testdata/hmmrank/../../modules/data/delete_me/gappa/gappa_taxonomy.tsv` (the `nf-core/modules` test-datasets `gappa_taxonomy.tsv` fixture, keyed by the same sequence ids) appended after the id, GTDB single-file style: `>id taxonomy;string`.
Content is otherwise byte-identical -- only the header lines changed.

Used to test deriving taxonomy from `--refseqfile` FASTA headers instead of a separate `--taxonomy` file, without needing a new reference alignment/tree pair -- the existing `PF14720_seed.ft.LGCAT.newick` tree still applies unchanged.

### `testdata/phyloplace_input.csv` - `reftreename` column

The sample sheet gained an optional `reftreename` column, naming the reference tree a row places onto so that rows sharing one can be summarised jointly as well as individually.
The existing rows already fell into two such groups, and the column just labels them: `PF14720wt` and `PF14720mafft` both place onto `PF14720_seed.ft.LGCAT.newick` with the same taxonomy file and are grouped as `PF14720`, while `nuc_hmmer` and `nuc_mafft` both place onto `cyanos_16s.newick` with no taxonomy and are grouped as `cyanos16s`.
Each group deliberately mixes the `hmmer` and `mafft` alignment methods, since placements made with different alignment methods still have to merge onto the same reference tree.

`PF14720wo` is left ungrouped on purpose, even though it places onto the same tree as the other two `PF14720` rows: it declares no taxonomy where they declare one, and a group has to agree on a taxonomy for a joint classification to mean anything.
So the file covers a group with a taxonomy, a group without one, and a row that opts out.

Nothing else in the file changed, and the column is ignored by pipeline versions that do not know about it.
Consumed by nf-core/phyloplace's `test_phyloplace_input` profile.

## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).

[^1]: From [stackoverflow](https://stackoverflow.com/a/60846265/11502856)
