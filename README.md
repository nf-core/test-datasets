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

## Datasets for nf-core/metatdenovo

### MetaEuk eukaryotic ORF-calling test data (RPL28 locus)

Created to support development of MetaEuk-based eukaryotic ORF calling ([issue #459](https://github.com/nf-core/metatdenovo/issues/459)), same-contig locus-level consolidation ([issue #463](https://github.com/nf-core/metatdenovo/issues/463)), and cross-contig protein-cluster-level consolidation ([issue #460](https://github.com/nf-core/metatdenovo/issues/460)).
The pipeline's existing test data has no eukaryotic genomic (intron-containing) content, so it can't exercise splice-aware gene calling or cross-caller/cross-contig consolidation.

Source: the *Saccharomyces cerevisiae* S288C RPL28 locus (60S ribosomal protein uL15), which has a single, well-annotated intron.
Downloaded from NCBI:

* Genomic region `NC_001139.9:310767-312127` (chromosome VII, includes the 511 bp intron plus flanking sequence)
* Spliced mRNA `NM_001180968.1` (exact CDS, no UTR)
* Protein `NP_011412.1` (149 aa)

Three read sets were simulated from these references with `wgsim -e 0.005 -r 0.0 -R 0.0` (sequencing error only, no true mutations/indels), 200 read pairs each:

* `rpl28_genomic` -- reads from the genomic (intron-containing) sequence.
* `rpl28_mrna` -- reads from the original spliced mRNA.
* `rpl28_mrna_recoded` -- reads from a synonymously recoded version of the same CDS (every codon substituted for its maximally-divergent synonymous alternative via a full genetic-code lookup table), 38.9% nucleotide-divergent from the original mRNA but translating to an identical 149 aa protein.

Assembly behaviour (MEGAHIT, `community.wave.seqera.io/library/megahit_pigz:87a590163e594224`), verified empirically before committing this data:

* `rpl28_genomic` + `rpl28_mrna` co-assemble into a single 1348 bp contig.
  MEGAHIT's graph merges the intron-skipping and intron-containing paths, since they share enough sequence identity/coverage support.
  This is the intended same-contig, splice-aware ORF-calling case for #459/#463: one contig, two valid gene models (spliced vs. unspliced) that a splice-aware caller like MetaEuk must resolve as one locus.
* `rpl28_genomic` + `rpl28_mrna_recoded` assemble into two separate contigs (1335 bp genomic, 449 bp recoded).
  The nucleotide divergence is large enough that MEGAHIT does not merge them.
* All three read sets combined assemble into exactly two contigs (1348 bp merged genomic+original-mRNA, 449 bp recoded-only).
  This is the intended cross-contig, protein-level consolidation case for #460: two contigs, two ORF calls, but one shared protein once translated.

`rpl28_reference_db.faa` (copy of `NP_011412.1`) is included as a minimal homology-search reference database for MetaEuk.

```
test_data/metaeuk/rpl28_genomic_R1.fastq.gz
test_data/metaeuk/rpl28_genomic_R2.fastq.gz
test_data/metaeuk/rpl28_mrna_R1.fastq.gz
test_data/metaeuk/rpl28_mrna_R2.fastq.gz
test_data/metaeuk/rpl28_mrna_recoded_R1.fastq.gz
test_data/metaeuk/rpl28_mrna_recoded_R2.fastq.gz
test_data/metaeuk/rpl28_reference_db.faa
```

### Single-row `diamond_dbs.csv` variant for CI-lightweight multi-caller testing

Added to support [issue #462](https://github.com/nf-core/metatdenovo/issues/462) (running multiple `--orf_caller` values in one execution).
The pipeline's `tests/megahit_prokka_transdecoder.nf.test` needs `--diamond_dbs` on to catch a real class of bug (a Nextflow `.join()` silently dropping duplicate left-side keys when multiple ORF callers share the same diamond db name) -- but that test already combines two ORF callers' own containers, and adding the full two-row `samplesheet/diamond_dbs.csv` (pulling both the `ncbi-refseq-test` and `gtdb-r220-test` databases/containers) pushed the combination over the CI runner's disk budget.

`diamond_dbs_single.csv` is the same `ncbi-refseq-test` row already used in `diamond_dbs.csv` (chosen over `gtdb-r220-test` since it also exercises the `parse_with_taxdump` code path), just as the only row -- still enough to reproduce the duplicate-key bug (two callers sharing one db name), at roughly half the diamond-side CI cost. No new database files -- points at the same `.dmnd`/`.dmp` files `diamond_dbs.csv` already references.

```
samplesheet/diamond_dbs_single.csv
```

[^1]: From [stackoverflow](https://stackoverflow.com/a/60846265/11502856)
