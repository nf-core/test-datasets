# test-datasets: `viralrecon`

This branch contains test data to be used for automated testing with the [nf-core/viralrecon](https://github.com/nf-core/viralrecon) pipeline.

## Content of this repository

### `samplesheet_test_sispa.csv`

Sample information sheet required to test the pipeline containing sample information and links to FastQ files stored in this repository. This sample sheet corresponds to Illumina SISPA data.

### `samplesheet_test_amplicon.csv`

Sample information sheet required to test the pipeline containing sample information and links to FastQ files stored in this repository. This sample sheet corresponds to Illumina amplicon data.

### `samplesheet_full_sispa.csv`

Sample information sheet required to test the pipeline containing sample information and links to original full FastQ files. This sample sheet corresponds to Illumina SISPA data.

### `samplesheet_full_amplicon.csv`

Sample information sheet required to test the pipeline containing sample information and links to original full FastQ files. This sample sheet corresponds to Illumina amplicon data.

### `genome/`

* `NC_045512.2.fasta`: Reference SARS-Cov2 fasta file.
* `NC_045512.2.gff`: Reference SARS-Cov2 gff file.

### `amplicon/nCoV-2019.artic.primers.fasta`

`fasta` file with the localization of the primers in the SARS-Cov-2 virus genome from an enrichment experiment using the Artic network amplicons. Retrieved from [Zenodo](https://doi.org/10.5281/zenodo.3735110).

### `amplicon/nCoV-2019.schemeMod.bed`

`bed` file with the localization of the primers in the SARS-Cov-2 virus genome from an enrichment experiment using the Artic network amplicons. Retrieved from [Zenodo](https://doi.org/10.5281/zenodo.3735110).

### `kraken2/kraken2_hs22.tar.gz`

Small host DB for kraken2 required to test the pipeline containing only human chr22. The commands used to generate the DB are:

```
# Download hs chr22
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz .
gunzip Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz

# Generate the DB
kraken2-build --db kraken2_hs22 --download-taxonomy
kraken2-build --db kraken2_hs22 --add-to-library Homo_sapiens.GRCh38.dna.chromosome.22.fa
kraken2-build --db kraken2_hs22 --build
```


### `fastq/illumina_sispa/`

| file                    | num_seqs | sum_len    | min_len | avg_len | max_len | file_size | Sequencer   | LibrarySource      |
|-------------------------|----------|------------|---------|---------|---------|-----------|-------------|--------------------|
| SRR11140744_R1.fastq.gz |   10,092 |  2,284,737 |     100 |   175.5 |     251 |      747K | PE Illumina | Metagenomics       |
| SRR11140744_R2.fastq.gz |   10,092 |  2,260,970 |     100 |   175.5 |     251 |      783K | PE Illumina | Metagenomics       |
| SRR11140746_R1.fastq.gz |    7,196 |  1,609,884 |     100 |   175.5 |     251 |      554K | PE Illumina | Metagenomics       |
| SRR11140746_R2.fastq.gz |    7,196 |  1,594,703 |     100 |   175.5 |     251 |      580K | PE Illumina | Metagenomics       |
| SRR11140748_R1.fastq.gz |    8,447 |  1,918,541 |     100 |   175.5 |     251 |      650K | PE Illumina | Metagenomics       |
| SRR11140748_R2.fastq.gz |    8,447 |  1,903,781 |     100 |   175.5 |     251 |      683K | PE Illumina | Metagenomics       |
| SRR11140750_R1.fastq.gz |      369 |     81,898 |     100 |   175.5 |     251 |       40K | PE Illumina | Metagenomics       |
| SRR11140750_R2.fastq.gz |      369 |     80,344 |     102 |   176.5 |     251 |       41K | PE Illumina | Metagenomics       |

> All FastQ files were sub-sampled to 0.02% of the original reads.

### `fastq/illumina_amplicon/`

| file                    | num_seqs | sum_len    | min_len | avg_len | max_len | file_size | Sequencer   | LibrarySource      |
|-------------------------|----------|------------|---------|---------|---------|-----------|-------------|--------------------|
| sample1_R1.fastq.gz     |   27,721 |  8,285,732 |      35 |     168 |     301 |        4M | PE Illumina | Metagenomics       |
| sample1_R2.fastq.gz     |   27,721 |  8,285,900 |      35 |     168 |     301 |        4M | PE Illumina | Metagenomics       |
| sample2_R1.fastq.gz     |   21,481 |  6,416,734 |      35 |     168 |     301 |        3M | PE Illumina | Metagenomics       |
| sample2_R2.fastq.gz     |   21,481 |  6,416,265 |      35 |     168 |     301 |        3M | PE Illumina | Metagenomics       |

> All FastQ files were sub-sampled to 0.02% of the original reads.

## Sampling procedure

Prepare a file `list.txt` with the following SRA accession numbers:

```
SRR11140744
SRR11140746
SRR11140748
SRR11140750
```

Download SRA using `parallel-fastq-dump` and `parallel`.

```bash
cat list.txt | parallel 'fasterq-dump {}'
```

Sub-sampling fastq files with a ratio of 0.02 using `seqkit`

```bash
parallel 'seqkit sample -p 0.02 -s 2020 {} | pigz > {.}.fastq.gz' ::: SRR*
```

The above tools are available on bioconda.

## Expected output

TBD.
