# test-datasets: `scnanoseq`

<!---
TODO: add link to scnanoseq pipeline when ready
-->
This branch contains test data to be used for automated testing with the nf-core/scnanoseq pipeline.

## Content of this repository

`reference/`: Sub-sampled genome reference files

`testdata/*.fastq.gz`: Minimal single-end gridion test dataset

`samplesheet/samplesheet.csv`: samplesheet file for minimal test dataset

`samplesheet/samplesheet_full.csv`: samplesheet file for full test dataset

## Minimal test dataset origin

*H. sapiens* single-end gridion single-cell transciptome data was obtained from:
>You, Y., Prawer, Y.D.J., De Paoli-Iseppi, R. et al. Identification of cell barcodes from long-read single-cell RNA-seq with BLAZE. Genome Biol 24, 66 (2023). https://doi.org/10.1186/s13059-023-02907-y

### Sampling information
| run_accession | experiment_alias | read_count | sample_title   |
|---------------|------------------|------------|----------------|
| ERR9958133    | ERX9501002       | 3423062    | Q20 Gridion    |
| ERR9958134    | ERX9501002       | 7521667    | LSK110 Gridion |

### Sampling procedure

Fastq files were subsetted to about 5000 reads using the below commands.
```
seqtk sample -s123 ERR9958133.fastq 5000 > ERR9958133_sub.fastq
seqtk sample -s123 ERR9958134.fastq 5000 > ERR9958134_sub.fastq

gzip ERR9958133_sub.fastq
gzip ERR9958133_sub.fastq
```

### Reference procedure

The test data in this repository is derived from *H. sapiens* data. Due to size of the reference files needed for this pipeline (fasta and gtf), they are beyond the scope of what is needed for testing these pipelines and bears unnecessary overhead for storage. To alleviate this burden, the reference files (which are GENCODE GRCh38 v40) have been filtered down to only contain information from chromosome 21.
