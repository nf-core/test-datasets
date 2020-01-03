# test-datasets: `slamseq`

This branch contains test data to be used for automated testing with the [nf-core/slamseq](https://github.com/nf-core/slamseq) pipeline.

## Content of this repository

`reference/`: Genome reference files (*Mus musculus* mm10 + [custom UTR annotation](https://github.com/AmeresLab/UTRannotation) subset to chr10)   
`testdata/` : FastQ files sub-sampled to 10,000 SR50 reads from chr10

## Dataset origin

*M. musculus* single-end 50 bp SLAM-seq dataset was obtained from:

Herzog VA, Reichholf B, Neumann T, Rescheneder P, Bhat P, Burkard TR, Wlotzka W, von Haeseler A, Zuber J & Ameres SL: Thiol-linked alkylation of RNA to assess expression dynamics. Nature Methods, 2017. 14(12), 1198–1204. [Pubmed](https://www.ncbi.nlm.nih.gov/pubmed/28945705) [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99978)

### Sample information

| GEO_ID	    | SRA_ID	    | SAMPLE_NAME	            |
|-------------|-------------|-------------------------|
| GSM2666816	| SRR5678869	| AN3-12 wt_mESC no s4U	  |
| GSM2666819	| SRR5678872	| AN3-12 wt_mESC 0h chase	|

## Expected output

To track and test the reproducibility of the pipeline with default parameters below are some of the expected outputs.
