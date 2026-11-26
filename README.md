# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

Test data to be used for automated testing with the nf-core pipelines

## Test data for Demultiplex

Data also used by modules lives on the
[`modules` branch](https://github.com/nf-core/test-datasets/tree/modules/data/genomics/homo_sapiens),
and the samplesheets on this branch point there when possible.

### Illumina Runfolders

#### MiSeq

`modules` branch: `data/genomics/homo_sapiens/illumina/bcl/flowcell_no_runparameters.tar.gz`
and `flowcell_samplesheet.csv`.

This MiSeq run has been trimmed to only keep the first tile in each cycle. To
demultiplex this data, you will need to use `--tiles s_1_1101` with `bcl2fasq`
or `--first-tile-only true` with `bclconvert`.

`modules` branch: `data/genomics/homo_sapiens/illumina/bcl/miseq_35147139_samplesheet.csv`
is the samplesheet for the full-size MiSeq run.

#### NextSeq 2000

`modules` branch: `data/genomics/homo_sapiens/illumina/bcl/nextseq2000_samplesheet.csv`.

#### NovaSeq 6000

This folder contains data from a NovaSeq 6000 run. The data has been trimmed to
only keep the first tile in each cycle. To demultiplex this data, you will need
to use `--tiles s_1_2101` with `bcl2fasq` or `--first-tile-only true` with
`bclconvert`.

### Element

`modules` branch: `data/genomics/homo_sapiens/element/sim-data/`
(`RunManifest.csv`, `out.sample_meta.csv`, `fastq.tar.gz`, `sim-data.tar.gz`,
`fqtk_sample_metadata_subset.tsv`).

### MGI

`modules` branch: `data/genomics/homo_sapiens/mgi/fc01_sample_sheet.csv`.
