# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

# test-datasets: `sammyseq`

This branch contains data to be used for automated testing with the [nf-core/sammyseq](https://github.com/daisymut/sammyseq) pipeline.

## Content of this repository

`testdata/CTRL004_S*_chr22only.fq.gz`: Human fibroblast single-end test data for pipeline sub-sampled to map on chr22 only

## Minimal test dataset origin

_H. sapiens_ fibroblast, 50bp single-end 3-fraction SAMMY-seq sequences was obtained from:

> Sebestyén, E., Marullo, F., Lucini, F. et al. SAMMY-seq reveals early alteration of heterochromatin and deregulation of bivalent genes in Hutchinson-Gilford Progeria Syndrome. Nat Commun 11, 6274 (2020). https://doi.org/10.1038/s41467-020-20048-9. [Pubmed](https://pubmed.ncbi.nlm.nih.gov/33293552/) [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118633)

### Sampling information

| GEO_sample | run_accession | read_count | SRA_experiment | sample_title         |
| ---------- | ------------- | ---------- | -------------- | -------------------- |
| GSM3335763 | SRR7610706    | 78683296   | SRX4475555     | CTRL004 SAMMY-seq S2 |
| GSM3335764 | SRR7610707    | 60438514   | SRX4475554     | CTRL004 SAMMY-seq S3 |
| GSM3335765 | SRR7610708    | 54864540   | SRX4475553     | CTRL004 SAMMY-seq S4 |
