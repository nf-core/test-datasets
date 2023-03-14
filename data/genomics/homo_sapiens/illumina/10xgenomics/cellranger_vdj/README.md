# cellranger vdj test data

This folder contains test data for the `cellranger vdj` module.

## Source

The test data are taken from the official 10X Genomics `cellranger vdj` tutorial:
- [Tutorial](https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/tutorial/tutorial-vdj)
- [Data](https://www.10xgenomics.com/resources/datasets/human-b-cells-from-a-healthy-donor-1-k-cells-2-standard-6-0-0)
- [Reference files](https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/tutorial/tutorial-vdj#download:~:text=https%3A//cf.10xgenomics.com/supp/cell%2Dvdj/refdata%2Dcellranger%2Dvdj%2DGRCh38%2Dalts%2Densembl%2D5.0.0.tar.gz)

Download the files as in the 10X tutorial:
```bash
# download a tarball of input FASTQ files
curl -LO https://cf.10xgenomics.com/samples/cell-vdj/6.0.0/sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex/sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex_fastqs.tar

# untar the FASTQs
tar -xf sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex_fastqs.tar

# download and untar the reference files
curl -O https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz
tar -xf refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz
```

Note that the tarball includes 5' gene expression FASTQs as well as B-cell sequencing FASTQs.
We only need the latter here.

## Subsampling

The original data are excessively large for nf-core testing purposes.
Decrease the file size by subsampling the reads. 
Note that `cellranger vdj` needs at least 10,000 reads to autodetect the library chemistry.
Data here were subsampled as follows:

```bash

# subsample B cell FASTQs in subfolder /sc5p_v2_hs_B_1k_b_fastqs/
for i in 1 2; do
  for j in 1 2; do
    zcat sc5p_v2_hs_B_1k_b_S1_L00${i}_R${j}_001.fastq.gz | head -n 40000 | gzip -c > subsampled_sc5p_v2_hs_B_1k_b_S1_L00${i}_R${j}_001.fastq.gz
  done
done

# rename subfolder to /subsampled_sc5p_v2_hs_B_1k_b_fastqs/
cd ..
mv sc5p_v2_hs_B_1k_b_fastqs/ subsampled_sc5p_v2_hs_B_1k_b_fastqs/
```

## Test data folder structure

```bash
cellranger_vdj/
├── README.md
├── refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0
│   ├── fasta
│   │   ├── regions.fa
│   │   └── supp_regions.fa
│   └── reference.json
└── subsampled_sc5p_v2_hs_B_1k_b_fastqs
    ├── subsampled_sc5p_v2_hs_B_1k_b_S1_L001_R1_001.fastq.gz
    ├── subsampled_sc5p_v2_hs_B_1k_b_S1_L001_R2_001.fastq.gz
    ├── subsampled_sc5p_v2_hs_B_1k_b_S1_L002_R1_001.fastq.gz
    └── subsampled_sc5p_v2_hs_B_1k_b_S1_L002_R2_001.fastq.gz
```
