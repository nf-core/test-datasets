# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

Test data to be used for automated testing with the nf-core pipelines

## Test data for Demultiplex

### Illumina Runfolders

#### MiSeq

This folder contains data from a MiSeq run. The data has been trimmed to only
keep the first tile in each cycle. To demultiplex this data, you will need to
use `--tiles s_1_1101` with `bcl2fasq` or `--first-tile-only true` with
`bclconvert`.
