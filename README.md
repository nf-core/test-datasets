# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

Test data to be used for automated testing with the nf-core pipelines

## Test data for Demultiplex

### Illumina Runfolders

#### MiSeq

This folder contains data from a MiSeq run. The data has been trimmed to only
keep the first tile in each cycle. To demultiplex this data, you will need to
use `--tiles s_1_1101` with `bcl2fasq` or `--first-tile-only true` with
`bclconvert`.

#### NovaSeq 6000

This folder contains data from a NovaSeq 6000 run. The data has been trimmed to
only keep the first tile in each cycle. To demultiplex this data, you will need
to use `--tiles s_1_2101` with `bcl2fasq` or `--first-tile-only true` with
`bclconvert`.

The MiSeq run and its `SampleSheet.csv` are also on the `modules` branch, as
`data/genomics/homo_sapiens/illumina/bcl/flowcell_no_runparameters.tar.gz` and
`flowcell_samplesheet.csv`.

### Samplesheets

`samplesheet/1.3.0/*_v2.csv` are copies of the matching samplesheets that point
to the identical files on the `modules` branch
(`data/genomics/homo_sapiens/element/sim-data/`,
`data/genomics/homo_sapiens/illumina/bcl/`). Data not on the `modules` branch
still comes from this branch, s3 or Hugging Face.
