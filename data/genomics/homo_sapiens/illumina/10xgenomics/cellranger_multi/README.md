# cellranger multi test data

This folder contains test data for the `cellranger multi` module.
Tests for `cellranger multi` will also refer to FASTQs from B cells in the [sister folder](../cellranger_vdj/README.md) for the `cellranger vdj` module.
This folder contains the corresponding 5' gene expression FASTQs for the B cell data.

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

# don't forget to download a transcriptome reference, e.g. from 10X Genomics
# or make a custom one with cellranger mkref
```

Note that the tarball includes 5' gene expression FASTQs as well as B-cell sequencing FASTQs.
We only need the *former* here.

Gene expression analyses in `cellranger multi` require a corresponding transcriptome reference.
Suitable ones can be downloaded [here](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest).
For tests, it's wisest to build a dummy reference against human chromosome 22 using the following files:
* [a FASTA reference](../../../genome/genome.fasta)
* [a corresponding GTF annotation](../../../genome/genome.gtf)

## Subsampling

The original data are excessively large for nf-core testing purposes.
Decrease the file size by subsampling the reads. 
Note that `cellranger` typically needs at least 10,000 reads to autodetect the library chemistry.
Since every FASTQ entry uses 4 lines, subsampling amounts to reading the first 40,000 lines of each FASTQ file.
Data here were subsampled as follows:

```bash

# subsample GEX FASTQs in subfolder /sc5p_v2_hs_B_1k_5gex_fastqs/
# some Unix systems might prefer gunzip in lieu of zcat
for i in 1 2; do
  for j in 1 2; do
    zcat sc5p_v2_hs_B_1k_5gex_S1_L00${i}_R${j}_001.fastq.gz | head -n 40000 | gzip -c > subsampled_sc5p_v2_hs_B_1k_5gex_S1_L00${i}_R${j}_001.fastq.gz
  done
done

# rename subfolder to /subsampled_sc5p_v2_hs_B_1k_5gex_fastqs/
cd ..
mv sc5p_v2_hs_B_1k_5gex_fastqs/ subsampled_sc5p_v2_hs_B_1k_5gex_fastqs/
```

While two lanes are provided in the original dataset,
one lane of data is sufficient for nf-core module tests.
We therefore discard lane 2.

## Test data folder structure

```bash
cellranger_multi/
├── README.md
└── subsampled_sc5p_v2_hs_B_1k_5gex_fastqs
    ├── subsampled_sc5p_v2_hs_B_1k_5gex_S1_L001_R1_001.fastq.gz
    └── subsampled_sc5p_v2_hs_B_1k_5gex_S1_L001_R2_001.fastq.gz
```
