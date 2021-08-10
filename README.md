# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

# test-datasets: `cutandrun`

Test data to be used for automated testing with the nf-core pipelines

This branch contains test data for the [nf-core/cutandrun](https://github.com/nf-core/cutandrun) pipeline.

## Documentation
nf-core/test-datasets comes with documentation in the `docs/` directory:

01. [Add a new  test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/ADD_NEW_DATA.md)
02. [Use an existing test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/USE_EXISTING_DATA.md)

## Support

For further information or help, don't hesitate to get in touch on the [Slack `#cutandrun` channel](https://nfcore.slack.com/channels/cutandrun) (you can join with [this invite](https://nf-co.re/join/slack)).

## Content of this repository

- `reference/`: Genome reference files used in testing
- `samplesheet/`: Sample sheets used for the `--input` param in testing
- `testdata/`: Experimental test data referenced in the samplesheets

## GSE145187 dataset origin

The test data available was taken from the original CUT&Tag protocol published in [nature](https://www.nature.com/articles/s41596-020-0373-x)/[pubmed](https://pubmed.ncbi.nlm.nih.gov/32913232/)

Kaya-Okur, H.S., Janssens, D.H., Henikoff, J.G. et al. Efficient low-cost chromatin profiling with CUT&Tag. Nat Protoc 15, 3264–3283 (2020). https://doi.org/10.1038/s41596-020-0373-x

The experiment was performed using two antibody targets for the histone marks h3k4me3 and h3k27me3. Each target had two biological replicates with an IgG control for each replicate.

The data can be found at [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145187)

### Folder Structure

- `GSE145187`: Contains the full test data from the GEO repository
- `GSE145187_10000`: Contains a random subset of 10,000 reads from the main dataset for testing