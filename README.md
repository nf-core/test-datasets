# test-datasets: `scflow`

This branch contains test data to be used for automated testing with the [nf-core/scflow](https://github.com/nf-core/scflow) pipeline.

## Content of this repository

`assets/ensembl_mappings.tsv`: A tsv file with biomart mappings between human ensembl_gene_id, gene_biotype, external_gene_name.   

`testdata/individual_*`: Directories containing matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz for individual samples.

`refs/Manifest.txt`: A tab-separated-variable file for the test dataset with two columns: key and filepath.  
`refs/SampleSheet.tsv`: A tab-separated-variable file with sample metadata for the test dataset.

## Minimal test dataset origin

Detailed information on how the minimal downsampled test data was generated can be found in [neurogenomics/scFlowExample](https://github.com/neurogenomics/scFlowExample) repo.
