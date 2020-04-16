# test-datasets: `viralrecon`

This branch contains test data to be used for automated testing with the [nf-core/viralrecon](https://github.com/nf-core/viralrecon) pipeline.

## Content of this repository

### `samplesheet/`

#### `samplesheet_test_sispa.csv`

Sample information sheet required to test the pipeline containing sample information and links to FastQ files stored in this repository. This sample sheet corresponds to Illumina SISPA data.

#### `samplesheet_test_amplicon.csv`

Sample information sheet required to test the pipeline containing sample information and links to FastQ files stored in this repository. This sample sheet corresponds to Illumina amplicon data.

#### `samplesheet_full_sispa.csv`

Sample information sheet required to test the pipeline containing sample information and links to original full FastQ files. This sample sheet corresponds to Illumina SISPA data.

#### `samplesheet_full_amplicon.csv`

Sample information sheet required to test the pipeline containing sample information and links to original full FastQ files. This sample sheet corresponds to Illumina amplicon data.

#### `samplesheet_test_sra.csv`
Sample information sheet required to test the pipeline containing sample information of one link to a original full FastQ files and two files that must be downloaded from SRA one single-end and one paired-end, respectively. This sample sheet corresponds to Illumina SISPA data.

### `genome/`

#### `kraken2/kraken2_hs22.tar.gz`

Small host database for `kraken2` containing only human chr22 required to test the pipeline. The commands used to generate the database are:

```
# Download human chr22
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz .
gunzip Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz

# Generate the database
kraken2-build --db kraken2_hs22 --download-taxonomy
kraken2-build --db kraken2_hs22 --add-to-library Homo_sapiens.GRCh38.dna.chromosome.22.fa
kraken2-build --db kraken2_hs22 --build
```

#### `NC_045512.2/`

* `GCF_009858895.2_ASM985889v3_genomic.<DOWNLOAD_DATE>.fna.gz`: SARS-CoV2 genome fasta file downloaded directly via [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz)
* `GCF_009858895.2_ASM985889v3_genomic.<DOWNLOAD_DATE>.gff.gz`: SARS-CoV2 genome GFF3 annotation file downloaded directly via [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz)
* `GCF_009858895.2_ASM985889v3_genomic.<DOWNLOAD_DATE>.gtf.gz`: SARS-CoV2 genome GTF2.2 annotation file downloaded directly via [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gtf.gz)
* `amplicon/`: ARTIC [V1](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V1), [V2](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V2) and [V3](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V3) primer schema files relative to the NC_045512.2 assembly. Files ending in `*.primer.fasta` were generated from the `.tsv` files in the repo.

#### `MN908947.3/`

* `GCA_009858895.3_ASM985889v3_genomic.200409.fna.gz`: SARS-CoV2 genome fasta file downloaded directly via [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/858/895/GCA_009858895.3_ASM985889v3/GCA_009858895.3_ASM985889v3_genomic.fna.gz)
* `GCA_009858895.3_ASM985889v3_genomic.200409.gff.gz`: SARS-CoV2 genome GFF3 annotation file downloaded directly via [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/858/895/GCA_009858895.3_ASM985889v3/GCA_009858895.3_ASM985889v3_genomic.gff.gz)
* `GCA_009858895.3_ASM985889v3_genomic.200409.gtf.gz`: SARS-CoV2 genome GTF2.2 annotation file downloaded directly via [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/858/895/GCA_009858895.3_ASM985889v3/GCA_009858895.3_ASM985889v3_genomic.gtf.gz)
* `amplicon/`: ARTIC [V1](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V1), [V2](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V2) and [V3](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V3) primer schema files relative to the NC_045512.2 assembly. Files ending in `*.primer.fasta` were generated from the `.tsv` files in the repo.

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
