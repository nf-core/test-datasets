# Cellranger references

This directory contains downloaded references for use in testing Cellranger modules and pipelines.

# References

## VDJ

For use with calls to `cellranger vdj` or `cellranger multi` calls that include V(D)J libraries.

### Source

The reference is found in the official 10X Genomics `cellranger vdj` tutorial:
- [Tutorial](https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/tutorial/tutorial-vdj)
- [Reference files](https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/tutorial/tutorial-vdj#download:~:text=https%3A//cf.10xgenomics.com/supp/cell%2Dvdj/refdata%2Dcellranger%2Dvdj%2DGRCh38%2Dalts%2Densembl%2D5.0.0.tar.gz)

```bash
# download and untar the reference files
curl -O https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz
tar -xf refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz
```
