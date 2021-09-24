# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

# test-datasets: `cutandrun`

This branch contains test data to be used for automated testing with the [nf-core/cutandrun](https://github.com/nf-core/cutandrun) nf-core pipeline. 

## Adding data to test-datasets

To add data to any branch on in nf-core/test-datasets follow the documentation below:

01. [Add a new test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/ADD_NEW_DATA.md)
02. [Use an existing test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/USE_EXISTING_DATA.md)

## Support

For further information or help, don't hesitate to get in touch on the [Slack `#cutandrun` channel](https://nfcore.slack.com/channels/cutandrun) (you can join with [this invite](https://nf-co.re/join/slack)).

## Content of this repository

- `reference/`: Human genome reference files subsampled to specific chromosomes
- `samplesheet/`: Sample sheets used for the `--input` param in testing
- `testdata/GSE145187`: Selected paried end samples
- `testdata/GSE145187_10000`: Selected paried end samples subsampled to 10,000 reads

## GSE145187 dataset origin

The test data available was taken from the original CUT&Tag protocol published in [nature](https://www.nature.com/articles/s41596-020-0373-x)/[pubmed](https://pubmed.ncbi.nlm.nih.gov/32913232/)

> Kaya-Okur, H.S., Janssens, D.H., Henikoff, J.G. et al. Efficient low-cost chromatin profiling with CUT&Tag. Nat Protoc 15, 3264–3283 (2020). https://doi.org/10.1038/s41596-020-0373-x

The experiment was performed using two antibody targets for the histone marks h3k4me3 and h3k27me3. Each target had two biological replicates with an IgG control for each replicate.

The data can be found at [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145187)

Specifically, the following datasets were used:

H3K27me3: 
  - SH_Hs_K27m3_NX_0918 as replicate 1: GEO accession: GSE145187, SRA entry: SRX8754646
  - SH_Hs_K27m3_Xpc_0107 as replicate 2: GEO accession: GSE145187, SRA entry: SRX7713678

H3K4me3:
  - SH_Hs_K4m3_NX_0918 as replicate 1: GEO accession: GSE145187, SRA entry: SRX7713692
  - SH_Hs_K4m3_Xpc_0107 as replicate 2: GEO accession: GSE145187, SRA entry: SRX7713696

IgG:
  - SH_Hs_IgG_1x_0924 as replicate 1:GEO accession: GSE145187, SRA entry: SRX8468909
  - SH_Hs_IgG_20181224 as replicate 2: GEO accession: GSM3680227, SRA entry: SRX5545346
