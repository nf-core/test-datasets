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

## Datasets for nf-core/sativa

### raxtax prefilter module test fixtures

Small synthetic fixtures for the `RAXTAX`/`RAXTAXFILTER` local module tests (the raxtax-based prefilter added ahead of the pipeline's EPA-ng placement stage).

`raxtax_test.fasta`: a self-classification dataset for `RAXTAX`'s own test. Five short (80 bp) random DNA sequences forming two genuinely-consistent clusters (`SeqA1`/`SeqA2` labeled cluster A, `SeqB1`/`SeqB2` labeled cluster B) plus `SeqMisl`, whose sequence content is close to cluster A (5 mutations from `SeqA1`) but is deliberately labeled with cluster B's taxonomy. Mutation distances were chosen so `SeqA1`'s nearest neighbour stays `SeqA2` (not `SeqMisl`), keeping the genuine cluster A members correctly self-consistent. Verified against the real `raxtax` binary (with `--skip-exact-matches`, as the pipeline runs it) to confirm the expected classification: `SeqMisl` calls back to cluster A at confidence 1.0, everything else confidently confirms its own declared cluster.

`raxtaxfilter_test.tax`, `raxtaxfilter_test.out`, `raxtaxfilter_test.fasta`: a hand-authored fixture for `RAXTAXFILTER`'s own test (`raxtaxfilter_test.out` is a synthetic `raxtax.out`-format file, not produced by the real binary). Covers a correctly-classified sequence (`SeqOK`), a sequence needing best-hit selection among multiple candidate lines where the winning hit agrees at genus but disagrees at species (`SeqSpeciesMismatch`), and a sequence disagreeing at every rank (`SeqSevereMismatch`). Used to test both the best-hit-selection logic and the single-check-point (not cumulative) semantics of `--filter-rank`: at `--filter-rank 1` (species) `SeqSpeciesMismatch` is flagged, at `--filter-rank 2` (genus) it isn't, while `SeqSevereMismatch` is flagged at both.

### GTDB archaeal 16S dataset

121 real 16S rRNA gene sequences extracted from [GTDB release 220](https://data.gtdb.ecogenomic.org/releases/release220/220.0/genomic_files_all/ssu_all_r220.fna.gz) (`ssu_all_r220.fna.gz`, filtered to `d__Archaea`), covering 10 named species hand-picked to form a nested taxonomic relationship ladder:

- `Haloferax volcanii` / `Haloferax prahovense` -- a genus pair
- `Halorubrum distributum` -- same family (Haloferacaceae) as the pair, different genus
- `Haloarcula pellucida` -- same order (Halobacteriales), different family
- `Halorutilus salinus` -- same class (Halobacteria), different order
- `Methanosarcina barkeri` -- same phylum (Halobacteriota), different class
- `Saccharolobus islandicus` / `Saccharolobus solfataricus` -- a second genus pair, different phylum entirely (Thermoproteota)
- `Methanocaldococcus jannaschii` -- different phylum (Methanobacteriota_A); the first archaeal genome ever sequenced
- `Pyrococcus furiosus` -- different phylum (Methanobacteriota_B)

GTDB's own taxonomy was kept as-is, including several genuine `raxtax` self-classification disagreements found while building the set (e.g. all 3 `Pyrococcus furiosus` sequences classify as `Methanocaldococcus jannaschii` at full confidence, despite being full-length, non-truncated sequences) -- left in undoctored as honest hard cases rather than resolved one way or the other, since the panel is small enough that this may just reflect a small comparison pool. On top of the real data, one deliberate swapped-label pair was added as a guaranteed positive control: `DupHaloA` (real `H. volcanii` content, tagged with `S. islandicus`'s lineage) and `DupSulfoA` (the reverse), maximally divergent so the mismatch is unambiguous.

Provided in three forms:

- `gtdb_archaea_16s_unaligned.fasta` + `gtdb_archaea_16s.tax`: today's two-file pipeline input (unaligned FASTA + tab-separated taxonomy; needs aligning before use).
- `gtdb_archaea_16s_aligned.fasta`: the same 121 sequences aligned with MAFFT (`mafft --thread 4`), 2541 columns.
- `gtdb_archaea_16s_combined_unaligned.fasta` / `gtdb_archaea_16s_combined_aligned.fasta`: the same sequences/alignment in GTDB's own single-file convention (`>id taxonomy` header, no separate taxonomy file) -- built ahead of planned pipeline support for that as an alternative to the two-file input.

Sequence IDs are kept as GTDB's native `ACCESSION~CONTIG` form. Used by the `test_gtdb` Nextflow profile (`conf/test_gtdb.config` in the pipeline repo) and its pipeline-level nf-tests, which exercise the raxtax prefilter on real, full-length sequence data rather than the small structural fixtures used elsewhere.

### Barrnap archaeal rRNA HMM database

`barrnap_arc.hmm`: [Barrnap](https://github.com/tseemann/barrnap)'s bundled archaeal rRNA HMM database, extracted unmodified from `/usr/local/lib/barrnap/db/arc.hmm` inside the `quay.io/biocontainers/barrnap:0.9--hdfd78af_4` container (`docker cp`, no rebuild). A multi-profile HMMER3 file containing four separately-named profiles (`16S_rRNA` [Rfam RF01959], `23S_rRNA`, `5S_rRNA` [Rfam RF00001], `5_8S_rRNA` [Rfam RF00002]).

Used as the default `--hmm` value (with `--hmm_name 16S_rRNA`) for the pipeline's unaligned-input support: unaligned sequences are aligned via `hmmalign` against a single profile fetched from this database with `hmmfetch`, before continuing through the rest of the pipeline as a normal alignment. Paired with the `gtdb_archaea_16s_unaligned.fasta` fixture above (same sequences as the pre-aligned GTDB dataset, letting the two be cross-checked against each other) in the `test_gtdb_unaligned` profile.

### CHECKALIGNED module test fixture

`checkaligned_test_unaligned.fasta`: three short synthetic sequences of different lengths (16/10/20 bp), used by `CHECKALIGNED`'s own test to confirm it correctly classifies non-uniform-length FASTA as unaligned (paired with the already-uniform-length `sativa_loo_test.fasta` fixture above for the aligned case).

### Gap filter module test fixtures

`gapfilter_test.fasta` + `gapfilter_test.tax`: four short synthetic sequences for `GAPFILTER`'s own test -- `SeqGood1` (fully populated, non-gap proportion 1.0), `SeqGood2` (non-gap proportion exactly 0.8), `SeqBad1` (0.3) and `SeqBad2` (0.0, all gaps). Covers both the pipeline's default `--min_nongap` threshold (0.3) and a stricter one (0.8), including the exactly-at-threshold boundary case (`>=`, not `>`).

### Embedded-taxonomy module test fixtures

`resolvetaxonomy_test_embedded.fasta`: two short synthetic sequences (`SeqA`, `SeqB`) with taxonomy embedded directly in their FASTA header, GTDB-style (`>id taxonomy;string`), for `RESOLVETAXONOMY`'s own test of deriving a taxonomy file from `--sequences` headers when no separate `--taxonomy` file is given.

`resolvetaxonomy_test.tax`: a separate taxonomy file for the same two sequence names, but with deliberately _different_ (and clearly distinguishable -- `Archaea` rather than `Bacteria`) lineages than what's embedded in `resolvetaxonomy_test_embedded.fasta`'s headers. Used to test that an explicitly provided `--taxonomy` file always wins over embedded header text, with a warning rather than silently using the wrong source.

## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).

[^1]: From [stackoverflow](https://stackoverflow.com/a/60846265/11502856)
