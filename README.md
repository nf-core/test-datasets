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

## Data description

Index of every file on this branch, grouped by the directory it lives in, with the [nf-core/viralmetagenome](https://github.com/nf-core/viralmetagenome) parameter it is passed through.

### `samplesheets/`

| File | Description | Passed via |
| --- | --- | --- |
| `samplesheet.csv` | Three paired-end samples (`SRR11140744`, `SRR11140748` and an `empty-SRR` entry) pointing at the `illumina/` FASTQs by raw GitHub URL. | `--input` (`test`, `test_fail_db`, `test_fail_mapped`, `test_minimal`, `test_umi` profiles) |
| `samplesheet_group.csv` | Same reads as above plus a `group` column that puts two samples in `group1` and one in `group2`, to exercise per-group merging. | `--input` (`test_group` profile) |
| `samplesheet_full.csv` | All four SRR accessions, used for the full-size run. | `--input` (`test_full` profile) |
| `metadata_test.tsv` | ENA run metadata (accessions, library layout, instrument, read counts, …) for the two real samples in `samplesheet.csv`, surfaced in the MultiQC report. | `--metadata` |
| `metadata_testfull.tsv` | The same ENA metadata table covering all four samples of `samplesheet_full.csv`. | `--metadata` (`test_full` profile) |
| `mapping_constraints.csv` | Three constraint rows (two fixed references, one `selection: true` row) whose `sequence`/`gff` columns point at `genomes/MN908947.3.*` and `db/sarscov2.fasta`. | `--mapping_constraints` (`test`, `test_umi` profiles) |
| `mapping_constraints_group.csv` | Same constraints keyed on group names instead of sample names, and without a `gff` column. | `--mapping_constraints` (`test_group` profile) |
| `mapping_constraints_fail.tsv` | A single adenovirus constraint that deliberately matches none of the SARS-CoV-2 samples, to check the pipeline fails gracefully. | `--mapping_constraints` (`test_fail_mapped` profile) |

### `illumina/`

| File(s) | Description | Passed via |
| --- | --- | --- |
| `SRR11140744_R{1,2}.fastq.gz` | ~10k paired reads subsampled from SARS-CoV-2 metagenomic run SRR11140744 (veroSTAT-1KO, Illumina MiSeq, PRJNA607948). | `fastq_1`/`fastq_2` of every samplesheet |
| `SRR11140746_R{1,2}.fastq.gz` | ~7k paired reads from SRR11140746, only used by the full-size run. | `fastq_1`/`fastq_2` of `samplesheet_full.csv` |
| `SRR11140748_R{1,2}.fastq.gz` | ~8k paired reads from SRR11140748 (vero76). | `fastq_1`/`fastq_2` of every samplesheet |
| `SRR11140750_R{1,2}.fastq.gz` | ~370 paired reads from SRR11140750 (swab), a deliberately shallow sample. | `fastq_1`/`fastq_2` of `samplesheet_full.csv` |
| `empty_{1,2}.fastq.gz` | Valid but empty gzipped FASTQs, so the pipeline is tested against a sample that yields no reads. | `fastq_1`/`fastq_2` of the `empty-SRR` row in `samplesheet.csv` |
| `HIV.contigs.fasta` | 15 HIV-1 SPAdes/MEGAHIT contigs used as the contig input when unit-testing the local `BLAST_FILTER` module. | nf-test input of `modules/local/blast_filter` |

### `genomes/`

| File | Description | Passed via |
| --- | --- | --- |
| `MN908947.3.fasta` | The SARS-CoV-2 Wuhan-Hu-1 reference genome. | `sequence` column of `mapping_constraints{,_group}.csv` |
| `MN908947.3.gff3` | Matching NCBI GFF3 annotation for Wuhan-Hu-1, used to annotate called variants. | `gff` column of `mapping_constraints.csv` |
| `hadv-a.fasta` | Human adenovirus A (`NC_001460.1`) — an intentionally unrelated genome. | `sequence` column of `mapping_constraints_fail.tsv` |
| `chr22_23800000-23980000.fa` | A 180 kb slice of human chromosome 22, the source sequence the `db/kraken2_hs22.tar.gz` host database was built from. | Not passed directly; kept for rebuilding `--host_k2_db` |
| `coverages_spades.idxstats` | A tiny samtools `idxstats` table giving per-contig coverage for SPAdes contigs. | nf-test input of the `fasta_contig_clust` subworkflow |

### `db/`

| File(s) | Description | Passed via |
| --- | --- | --- |
| `sarscov2.fasta` | 135 SARS-CoV-2 genomes acting as the pool of candidate references for the contigs. | `--reference_pool` |
| `ebov-zaire.fasta` | 20 Zaire ebolavirus genomes — a reference pool that deliberately matches nothing in the SARS-CoV-2 test reads. | `--reference_pool` (`test_fail_db` profile) |
| `C-RVDBv29.0.fasta.subset.gz` | A 1,000-sequence subset of the clustered nucleotide Reference Viral DataBase v29.0, mirroring the pipeline's default reference pool. Not currently wired into any test profile. | `--reference_pool` |
| `kraken2_hs22.tar.gz` | Minimal Kraken2 database built from the human chr22 slice in `genomes/`, for host/contamination removal. | `--host_k2_db` |
| `kraken2_bracken.tar.gz` | Minimal Kraken2 database with the Bracken `*.kmer_distrib` files included, for read classification and abundance re-estimation. | `--kraken2_db` |
| `kaiju.tar.gz` | Minimal Kaiju database (`proteins.fmi` plus NCBI `names.dmp`/`nodes.dmp`) for protein-level read classification. | `--kaiju_db` |
| `virosaurus90_vertebrate-20200330.subset.fas.gz` | 1,000 sequences from the Virosaurus 90 vertebrate release, used to annotate the consensus constructs. | `--annotation_db` |
| `U-RVDBv29.0-prot_clustered.subset.fasta.xz` | 1,000 viral protein sequences from the clustered protein RVDB v29.0, handed to Prokka as a custom `--protein` database. | `--prokka_db` (`test_full` profile) |
| `checkv_minimal_db.tar` | A stripped-down CheckV database (`genome_db` + `hmm_db`) for consensus quality control. | `--checkv_db` (`test_full` profile) |
| `HIV.fasta` | 26 HIV-1 reference sequences used to build the BLAST database in the local `BLAST_FILTER` unit test. | nf-test input of `modules/local/blast_filter` |
| `sarscov2-uk-ncbi-virus.fasta.gz` | 74 UK SARS-CoV-2 genomes downloaded from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/), demonstrating a custom annotation database built outside Virosaurus ([viralmetagenome#318](https://github.com/nf-core/viralmetagenome/pull/318)). | `--annotation_db` |
| `sarscov2-uk-ncbi-virus.metadata.csv` | The metadata table exported alongside those 74 genomes (accession, isolate, collection date, geo-location, Pangolin lineage, host, …). The FASTA headers carry no embedded `key="value"` annotations, so the annotations are read from this table instead — it must be supplied together with the FASTA. | `--annotation_metadata` (alongside `--annotation_db`) |

### `local_modules/`

| File | Description | Passed via |
| --- | --- | --- |
| `blast_filter/HIV.blast.out.txt` | A pre-computed BLAST tabular result for `illumina/HIV.contigs.fasta` against `db/HIV.fasta`, so the module can be tested without running BLAST. | nf-test input of `modules/local/blast_filter` |
| `blast_filter/HIV.blacklist.txt` | Four accessions to exclude from reference selection, plus one deliberately misspelled entry to check the module tolerates it. | `--blacklist` / nf-test input of `modules/local/blast_filter` |

> If you add new files to any of these directories, please extend the matching table in the same PR.

## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).

[^1]: From [stackoverflow](https://stackoverflow.com/a/60846265/11502856)
