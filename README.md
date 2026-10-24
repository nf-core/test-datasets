# test-datasets: `scdownstream`

This branch contains test data to be used for automated testing with the [nf-core/scdownstream](https://github.com/nf-core/scdownstream) pipeline.

## Content of this repository

- `pbmc/`: Biological count inputs and shared PCA, neighbour graphs and Leiden clusters.
- `samplesheet.csv` and `samplesheet_single.csv`: Pipeline input samplesheets using the upstream `nf-core/test-datasets` `scdownstream` branch.
- `extension_base/`: Matching compact AnnData, scVI and Symphony references with raw-count provenance.
- `volcanoplot/`: Synthetic differential-expression tables for plotting regressions and their generator.
- `anndata-variations/`, `intermediates/`, `misc/` and `singleR/`: Compact structural fixtures, shared derivatives, format inputs and annotation references.
- `generation/`, `fixtures.json` and `requirements.txt`: Generation recipes, fixture contracts and required Python packages.

Small text inputs, including custom decoupler networks, sample identifier validation sheets and interesting-gene lists, are generated inline by the pipeline's nf-test tests.

## Fixture generation and provenance

The biological panel retains 384 cells per source and 3,000 shared genes. Raw inputs also retain genuine background droplets for EmptyDrops and SCAR. Counts are subsetted, never inflated. Individual sources preserve source identities; combined inputs expose explicitly synthetic donor, condition and label fields for test design.

All 72 generated fixture files are below 1,000,000 bytes. H5AD arrays use gzip compression; identical `X` and `raw_counts` matrices share HDF5 storage. The two root samplesheets are additional metadata files. The previous fixture inventory had 91 files; the published inventory has 74 including those samplesheets.

Run generation into a new external directory, with Apptainer and the pinned pipeline checkout available:

```bash
python -m pip install -r requirements.txt
python -m generation build --stage all \
    --out /path/to/new-candidate --cache /path/to/source-cache \
    --pipeline /path/to/scdownstream
python -m generation validate \
    --dataset /path/to/new-candidate --pipeline /path/to/scdownstream \
    --report /path/to/validation.json
python -m generation audit \
    --dataset /path/to/new-candidate --pipeline /path/to/scdownstream \
    --report /path/to/consumer-audit.json
```

Generation writes stage reports to the external output directory; reports are not committed to this repository. The compatibility contract pins the baseline pipeline revision. Synthetic volcano tables retain plotting edge cases that biological differential-expression outputs do not reliably produce.

Synthetic fixtures and generation scripts use the repository's MIT licence. Biological derivatives retain the source provenance below.

## Pipeline migration

The companion pipeline branch is `test/rebuilt-test-datasets`. Dataset and reference URLs target the upstream `nf-core/test-datasets` `scdownstream` branch, where these files will be available after merge. The dataset changes must merge before the companion pipeline changes. `pbmc/combined_filtered_pca.h5ad` is the canonical labelled multi-source input; `extension_base/` is the matching reference family. Structural inputs retain missing-raw-count rejection coverage without another learned reference family. Previously published commit-pinned URLs remain available.

Validation covered 420 enabled cases: 418 passed across full category runs and focused no-update verification. All ten pipeline cases passed with snapshot updates disabled. Two live CyteType cases remain blocked by hosted-service authentication; they are not passing annotation coverage.

## Test dataset origin

The `SRR28679756` to `SRR28679759` input names identify runs from [PRJNA1100403](https://www.ebi.ac.uk/ena/browser/view/PRJNA1100403), also registered as [GSE263932](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263932): "Selective CAR-T cell mediated B cell depletion suppresses interferon signature in SLE [Bcells]". The source study describes isolated B cells. The historical `pbmc/` directory name does not establish a PBMC cell-type composition.

ENA records `SRR28679756` as patient 7 after treatment, `SRR28679757` and `SRR28679758` as patient 7 before treatment, and `SRR28679759` as patient 6 after treatment. Runs `SRR28679757` and `SRR28679758` share experiment `SRX24247853` and sample `SAMN40969802`; they are not independent biological donors. These accession records were checked on 2026-10-23. They establish the run identities, not the complete processing history of the existing fixture matrices.

Labels and conditions generated for tests are synthetic unless explicitly documented otherwise. Synthetic donor groups do not establish biological replication.

### Test dataset pre-processing

The pinned fixture blobs establish the source bytes. Their original processing and downsampling history are not independently verified; this rebuild does not claim a fresh `nf-core/scrnaseq` run.
