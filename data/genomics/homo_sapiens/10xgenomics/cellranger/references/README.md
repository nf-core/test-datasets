# Cellranger references

This directory contains downloaded references for use in testing Cellranger modules and pipelines.

## References

### GEX

This repository does not store references for `cellranger count` and related single cell gene expression analyses.
Instead, choose one of two options:
* either use `cellranger mkref` to construct one from an existing reference sequence, such as human [chromosome 21](../../../../genome/chr21/sequence/genome.fasta) or [chromosome 22](../../../../genome/genome.fasta); or
* [download](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build) a full reference from 10X Genomics

### VDJ

For use with calls to `cellranger vdj` or `cellranger multi` calls that include V(D)J libraries.
The option to [download a VDJ reference from IMGT](https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/advanced/references#imgt) is _not supported_ in nf-core but will otherwise work with a correctly configured `cellranger` container or installation.

#### Source

The reference is found in the official 10X Genomics `cellranger vdj` tutorial:
- [Tutorial](https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/tutorial/tutorial-vdj)
- [Reference files](https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/tutorial/tutorial-vdj#download:~:text=https%3A//cf.10xgenomics.com/supp/cell%2Dvdj/refdata%2Dcellranger%2Dvdj%2DGRCh38%2Dalts%2Densembl%2D5.0.0.tar.gz)

```bash
# download and untar the reference files
curl -O https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz
tar -xf refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz
```
