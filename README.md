# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the nf-core pipelines

> ⚠️ **Do not merge your test data to `master`! Each pipeline has a dedicated branch (and a special one for modules)**

## Introduction

nf-core is a collection of high quality Nextflow pipelines. This repository contains various files for CI and unit testing of nf-core pipelines and infrastructure.

The principle for nf-core test data is as small as possible, as large as necessary. Please see the [guidelines](https://nf-co.re/docs/contributing/test_data_guidelines) for more detailed information. Always ask for guidance on the [nf-core slack](https://nf-co.re/join) before adding new test data.

## Test datasets for nf-core/diseasemodulediscovery

This branch holds the test data for the [nf-core/diseasemodulediscovery](https://github.com/nf-core/diseasemodulediscovery) pipeline. The pipeline takes a set of *seed* genes/proteins together with a protein–protein interaction (PPI) *network* and infers a disease module around the seeds. The datasets below are organised accordingly into seed files (one identifier per line) and network files (comma-separated edge lists, one edge per line).

### `small/`

Minimal seed and network files for fast CI runs. They are adapted from the example data of the [DIAMOnD](https://github.com/dinaghiassian/DIAMOnD/tree/master/Example) repository (MIT license).

**Networks** (comma-separated edge lists of Entrez gene IDs):

| File | Description |
| --- | --- |
| `entrez_ppi.csv` | PPI network derived from the DIAMOnD example. |
| `entrez_ppi_small.csv` | Smaller variant of the above for even faster runs. |

**Seed files** — the same seed set is provided in several identifier formats so the pipeline's ID mapping can be tested against each:

| File | Identifier type |
| --- | --- |
| `entrez_seeds_1.csv` | Entrez gene ID |
| `entrez_seeds_2.csv` | Entrez gene ID |
| `entrez_single_seed.csv` | Entrez gene ID (edge case: single seed) |
| `ensembl_seeds_2.csv` | Ensembl gene ID |
| `symbol_seeds_2.csv` | HGNC gene symbol |
| `uniprot_seeds.csv` | UniProt accession |

#### `small/permuted_networks/`

Precomputed permuted (degree-preserving re-wired) versions of the `small` networks, stored in [graph-tool](https://graph-tool.skewed.de/) binary `.gt` format. They let permutation-based significance testing run without re-generating permutations on every CI run. Three permutations (`perm_0`, `perm_1`, `perm_2`) are provided for each network under `entrez_ppi/` and `entrez_ppi_small/`.

### `mmusculus/`

The `small` data converted to mouse (*Mus musculus*) so the pipeline can be tested with a non-human organism. Identifiers are Ensembl mouse gene IDs (`ENSMUSG…`).

| File | Description |
| --- | --- |
| `ENSMUSG_network.csv` | Mouse PPI network (comma-separated edge list). |
| `ENSMUSG_seeds.txt` | Mouse seed set. |

### `issue_149/`

Data created to reproduce and test the fix for [nf-core/diseasemodulediscovery#149](https://github.com/nf-core/diseasemodulediscovery/issues/149).

| File | Description |
| --- | --- |
| `string.human_physical_links_v12_0_min700.Entrez.csv` | Human STRING physical-interaction network (v12.0, combined score ≥ 700), mapped to Entrez gene IDs. |
| `mondo.0009025.tsv` | Single-seed input, named after its MONDO disease ID. |

### `full/`

Input for the full-size test profile that exercises the pipeline with a large seed set.

| File | Description |
| --- | --- |
| `LUAD.tsv` | Lung adenocarcinoma (LUAD) associated gene symbols. |

LUAD is one of the complex diseases evaluated in the pipeline publication; see Kersting *et al.*, *Inferring and evaluating network medicine-based disease modules with Nextflow*, Bioinformatics (2026), [doi:10.1093/bioinformatics/btag223](https://doi.org/10.1093/bioinformatics/btag223).


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
