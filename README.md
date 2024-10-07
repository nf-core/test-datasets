# test-datasets: `genomeqc`

This branch contains test data to be used for automated testing with the [nf-core/genomeqc](https://github.com/nf-core/genomeqc) pipeline.

Currently this is in a temporary repo:
[Eco-Flow/genomeqc](https://github.com/Eco-Flow/genomeqc).

## Content of this repository

- `testdata/` : subsampled genome fasta files and annotations (gff).
- `samplesheet/input_bacteria.csv`: Minimal test dataset, whole mycoplasma genomes. Just links to input Refseq IDs (which downloads the data in first process).
- `samplesheet/input_myco_tiny.csv`: Minimal test dataset (actual data in this repo), partial mycoplasma genomes and their annotations, taking only the first 300 line of the GFF (or nearest to include whole genes), for two species only (Mycoplasmoides_fastidiosum, Mycoplasmoides_genitalium).

## Test datasets origin

NCBI:

```bash
Mycoplasmoides_fastidiosum,GCF_024498275.1
Mycoplasmoides_gallisepticum,GCF_017654545.1
Mycoplasmoides_pneumoniae,GCF_000733995.1
Mycoplasmoides_genitalium,GCF_000027325.1
```

