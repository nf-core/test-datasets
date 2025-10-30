# test-datasets: `viralrecon`

This branch contains test data to be used for automated testing with the [nf-core/viralrecon](https://github.com/nf-core/viralrecon) pipeline.

## Content of this repository

### `samplesheet/`

#### `samplesheet_test_nanopore.csv`

Sample information sheet required to test the pipeline containing sample names and barcodes for MinION data hosted in this repository. For testing purposes, some barcodes have been appended to this samplesheet that may not necessarily have associated data. The raw data associated with this run can be found in [`nanopore/minion`](nanopore/minion), and has been sub-setted to include a maximum of 3 `fast5`/`fastq`files per barcode.

#### `samplesheet_test_sra.csv`

Sample information sheet required to test the pipeline containing sample information of one link to a original full FastQ files and two files that must be downloaded from SRA one single-end and one paired-end, respectively. This sample sheet corresponds to Illumina SISPA data.

#### `samplesheet_test_illumina_sispa.csv`

Sample information sheet required to test the pipeline containing sample information and links to FastQ files stored in this repository. This sample sheet corresponds to Illumina SISPA data.

#### `samplesheet_test_illumina_amplicon.csv`

Sample information sheet required to test the pipeline containing sample information and links to FastQ files stored in this repository. This sample sheet corresponds to Illumina amplicon data.

#### `samplesheet_full_illumina_sispa.csv`

Sample information sheet required to test the pipeline containing sample information and links to original full FastQ files. This sample sheet corresponds to Illumina SISPA data.

#### `samplesheet_full_illumina_amplicon.csv`

Sample information sheet required to test the pipeline containing sample information and links to original full FastQ files. This sample sheet corresponds to Illumina amplicon data.

#### `samplesheet_test_EV.csv`

Sample information sheet required to test the pipeline containing sample information and links to FastQ files stored in this repository. This sample sheet corresponds to subsampled Illumina metagenomics enterovirus data.

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

#### `blastdb/ev_test_blastdb.tar.gz`

Small custom blast database with containing 8 enterovirus genomes with taxid mapping required to test the pipeline. The steps used to generate the database are:

Prepare a fasta file `genomes.fasta` with genomes from the accession numbers below.

Prepare a file `taxid_map.txt` with one accession number and corresponding tab separated tax id per line:

```
U22521.1        39054
AY421764.1      86107
KC507895.1      31704
AF162711.1      41846
AF114383.1      12074
AF538841.1      12080
NC_038308.1     42789
JX393302.1      2749421
```

Create custom blast database using `makeblastdb`

```bash
makeblastdb -in genomes.fasta -dbtype nucl -parse_seqids -out ev_test_blastdb -taxid_map taxid_map.txt
```

Manually download the files `taxdb.bt[id]` from NCBI, through [this link](https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz/) and add to the created blast database folder

Compress using `tar`

```bash
tar -czvf ev_test_blastdb.tar.gz ev_test_blastdb/
```

#### `NC_045512.2/`

-   `GCF_009858895.2_ASM985889v3_genomic.<DOWNLOAD_DATE>.fna.gz`: SARS-CoV2 genome fasta file downloaded directly via [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz)
-   `GCF_009858895.2_ASM985889v3_genomic.<DOWNLOAD_DATE>.gff.gz`: SARS-CoV2 genome GFF3 annotation file downloaded directly via [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz)
-   `GCF_009858895.2_ASM985889v3_genomic.<DOWNLOAD_DATE>.gtf.gz`: SARS-CoV2 genome GTF2.2 annotation file downloaded directly via [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gtf.gz)
-   `amplicon/`: ARTIC [V1](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V1), [V2](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V2) and [V3](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V3) primer schema files relative to the NC_045512.2 assembly. Files ending in `*.primer.fasta` were generated from the `.tsv` files in the repo.

#### `MN908947.3/`

-   `GCA_009858895.3_ASM985889v3_genomic.<DOWNLOAD_DATE>.fna.gz`: SARS-CoV2 genome fasta file downloaded directly via [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/858/895/GCA_009858895.3_ASM985889v3/GCA_009858895.3_ASM985889v3_genomic.fna.gz)
-   `GCA_009858895.3_ASM985889v3_genomic.<DOWNLOAD_DATE>.gff.gz`: SARS-CoV2 genome GFF3 annotation file downloaded directly via [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/858/895/GCA_009858895.3_ASM985889v3/GCA_009858895.3_ASM985889v3_genomic.gff.gz)
-   `GCA_009858895.3_ASM985889v3_genomic.<DOWNLOAD_DATE>.gtf.gz`: SARS-CoV2 genome GTF2.2 annotation file downloaded directly via [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/858/895/GCA_009858895.3_ASM985889v3/GCA_009858895.3_ASM985889v3_genomic.gtf.gz)
-   `amplicon/`: ARTIC [V1](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V1), [V2](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V2) and [V3](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V3) primer schema files relative to the MN908947.3 assembly. Files ending in `*.primer.fasta` were generated from the `.tsv` files in the repo.
-   `nextclade_sars-cov-2_MN908947_2024-10-17--16_48_48Z.tar.gz`: A set of input data files required for Nextclade to run an analysis on SARS-CoV2. Previous format did not require `pathogen.json` file but from v3+ it is required. File was created with `nextclade dataset get -n sars-cov-2 --tag 2024-10-17--16-48-48Z`.

#### `NC_063383.1`

-   `GCF_014621545.1_ASM1462154v1_genomic.<DOWNLOAD_DATE>.fna.gz`: Monkeypox genome fasta file downloaded directly via [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/621/545/GCF_014621545.1_ASM1462154v1/GCF_014621545.1_ASM1462154v1_genomic.fna.gz)
-   `GCF_014621545.1_ASM1462154v1_genomic.<DOWNLOAD_DATE>.gff.gz`: Monkeypox genome GFF3 annotation file downloaded directly via [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/621/545/GCF_014621545.1_ASM1462154v1/GCF_014621545.1_ASM1462154v1_genomic.gff.gz)
-   `nextclade_hMPXV_NC_063383.1_2024-08-27--21-28-04Z.tar.gz`: A set of input data files required for Nextclade to run an analysis on MPOX. Previous format did not require `pathogen.json` file but from v3+ it is required. File was created with `nextclade dataset get -n MPXV --tag 2024-08-27--21-28-04Z`.

#### `ON563414.3`

-   `GCA_023516015.3_ASM2351601v1_genomic.<DOWNLOAD_DATE>.fna.gz`: Monkeypox genome fasta file downloaded directly via [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/023/516/015/GCA_023516015.3_ASM2351601v1/GCA_023516015.3_ASM2351601v1_genomic.fna.gz)
-   `GCA_023516015.3_ASM2351601v1_genomic.<DOWNLOAD_DATE>.gff.gz`: Monkeypox genome GFF3 annotation file downloaded directly via [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/023/516/015/GCA_023516015.3_ASM2351601v1/GCA_023516015.3_ASM2351601v1_genomic.gff.gz)

#### `MT903344.1`

-   `GCA_014621585.1_ASM1462158v1_genomic.<DOWNLOAD_DATE>.fna.gz`: Monkeypox genome fasta file downloaded directly via [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/621/585/GCA_014621585.1_ASM1462158v1/GCA_014621585.1_ASM1462158v1_genomic.fna.gz)
-   `GCA_014621585.1_ASM1462158v1_genomic.<DOWNLOAD_DATE>.gff.gz`: Monkeypox genome GFF3 annotation file downloaded directly via [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/621/585/GCA_014621585.1_ASM1462158v1/GCA_014621585.1_ASM1462158v1_genomic.gff.gz)

#### `NC_002058.3`

-   `NC_002058.3.fasta`: Enterovirus C: Human Poliovirus 1 genome fasta file downloaded directly via [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/861/165/GCF_000861165.1_ViralProj15288/GCF_000861165.1_ViralProj15288_genomic.fna.gz)
-   `NC_002058.3.fasta.gff`: Enterovirus C: Human Poliovirus 1 genome GFF3 annotation file downloaded directly via [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/861/165/GCF_000861165.1_ViralProj15288/GCF_000861165.1_ViralProj15288_genomic.gff.gz)

### `fastq/illumina_sispa/`

| file                    | num_seqs | sum_len   | min_len | avg_len | max_len | file_size | Sequencer   | LibrarySource |
| ----------------------- | -------- | --------- | ------- | ------- | ------- | --------- | ----------- | ------------- |
| SRR11140744_R1.fastq.gz | 10,092   | 2,284,737 | 100     | 175.5   | 251     | 747K      | PE Illumina | Metagenomics  |
| SRR11140744_R2.fastq.gz | 10,092   | 2,260,970 | 100     | 175.5   | 251     | 783K      | PE Illumina | Metagenomics  |
| SRR11140746_R1.fastq.gz | 7,196    | 1,609,884 | 100     | 175.5   | 251     | 554K      | PE Illumina | Metagenomics  |
| SRR11140746_R2.fastq.gz | 7,196    | 1,594,703 | 100     | 175.5   | 251     | 580K      | PE Illumina | Metagenomics  |
| SRR11140748_R1.fastq.gz | 8,447    | 1,918,541 | 100     | 175.5   | 251     | 650K      | PE Illumina | Metagenomics  |
| SRR11140748_R2.fastq.gz | 8,447    | 1,903,781 | 100     | 175.5   | 251     | 683K      | PE Illumina | Metagenomics  |
| SRR11140750_R1.fastq.gz | 369      | 81,898    | 100     | 175.5   | 251     | 40K       | PE Illumina | Metagenomics  |
| SRR11140750_R2.fastq.gz | 369      | 80,344    | 102     | 176.5   | 251     | 41K       | PE Illumina | Metagenomics  |

> All FastQ files were sub-sampled to 0.02% of the original reads.

### `fastq/illumina_amplicon/`

| file                | num_seqs | sum_len   | min_len | avg_len | max_len | file_size | Sequencer   | LibrarySource |
| ------------------- | -------- | --------- | ------- | ------- | ------- | --------- | ----------- | ------------- |
| sample1_R1.fastq.gz | 27,721   | 8,285,732 | 35      | 168     | 301     | 4M        | PE Illumina | Metagenomics  |
| sample1_R2.fastq.gz | 27,721   | 8,285,900 | 35      | 168     | 301     | 4M        | PE Illumina | Metagenomics  |
| sample2_R1.fastq.gz | 21,481   | 6,416,734 | 35      | 168     | 301     | 3M        | PE Illumina | Metagenomics  |
| sample2_R2.fastq.gz | 21,481   | 6,416,265 | 35      | 168     | 301     | 3M        | PE Illumina | Metagenomics  |

> All FastQ files were sub-sampled to 0.02% of the original reads.

### `fastq/illumina_enterovirus/`

| file                   | num_seqs | sum_len   | min_len | avg_len | max_len | file_size | Sequencer   | LibrarySource |
| ---------------------- | -------- | --------- | ------- | ------- | ------- | --------- | ----------- | ------------- |
| SRR13266665_1.fastq.gz | 10,000   | 1,342,331 | 35      | 134.23  | 151     | 508K      | PE Illumina | Metagenomics  |
| SRR13266665_2.fastq.gz | 10,000   | 1,321,544 | 35      | 132.15  | 151     | 508K      | PE Illumina | Metagenomics  |

> All FastQ files were sub-sampled down to 10,000 reads.

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

For enterovirus:
Sub-sampling fastq files to 10,000 reads using `seqtk`

```bash
seqtk sample -s100 SRR13266665_1.fastq 10000 > sub_SRR13266665_1.fastq
seqtk sample -s100 SRR13266665_2.fastq 10000 > sub_SRR13266665_2.fastq
```

The above tools are available on bioconda.

## Expected output

TBD.
