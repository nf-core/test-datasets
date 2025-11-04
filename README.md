# test-datasets: `viralrecon`

This branch contains test data to be used for automated testing with the [nf-core/viralrecon](https://github.com/nf-core/viralrecon) pipeline.

## Content of this repository

### `samplesheet/`

#### `samplesheet_test_nanopore.csv`

Sample information sheet required to test the pipeline containing sample names and barcodes for MinION data hosted in this repository. For testing purposes, some barcodes have been appended to this samplesheet that may not necessarily have associated data. The raw data associated with this run can be found in [`nanopore/minion`](nanopore/minion), and has been sub-setted to include a maximum of 3 `fast5`/`fastq`files per barcode.

#### `samplesheet_test_sra.csv`

Sample information sheet required to test the pipeline containing sample information of one link to a original full FastQ files and two files that must be downloaded from SRA one single-end and one paired-end, respectively.

This sample sheet corresponds to SARS-CoV-2 Illumina SISPA data.

#### `samplesheet_test_illumina_sispa.csv`

Sample information sheet required to test the pipeline containing sample information and links to FastQ files stored in this repository.

This sample sheet corresponds to SARS-CoV-2 Illumina SISPA data.

#### `samplesheet_test_illumina_amplicon.csv`

Sample information sheet required to test the pipeline containing sample information and links to FastQ files stored in this repository.

This sample sheet corresponds to SARS-CoV-2 Illumina amplicon primer enrichment data.

#### `samplesheet_full_illumina_sispa.csv`

Sample information sheet required to test the pipeline containing sample information and links to original full FastQ files.

This sample sheet corresponds to SARS-CoV-2 Illumina SISPA data.

#### `samplesheet_full_illumina_amplicon.csv`

Sample information sheet required to test the pipeline containing sample information and links to original full FastQ files.

This sample sheet corresponds to SARS-CoV-2 Illumina amplicon primer enrichment data.

#### `samplesheet_full_illumina_fragmented.csv`

Sample information sheet required to test the pipeline containing sample information and links to FastQ files stored in this repository.

This sample sheet corresponds to Crimea Congo data.

#### `v3.0/samplesheet_test_hiv.csv`

Sample information sheet required to test the pipeline containing sample information and links to FastQ files stored in this repository.

This sample sheet corresponds to HIV Illumina amplicon primer enrichment data from different SRA experiments. For test purposes these will be trated as non amplicon data.

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


#### `crimea_congo`

Crimea Congo is a fragmented genome with three fragments. S, M and L based on the fragment size.

-   `crimea_congo.fasta.gz`: Crimea Congo fasta genome containing S, M and L fragments: KY484036.1, KY484035.1, KY484034.1
-   `crimea_congo.gff.gz`: Crimea congo genome GFF3 annotation file containing annotation for S, M and L fragments: KY484036.1, KY484035.1, KY484034.1

#### `NC_001802.1`

This reference was chosen based on [Nextclade's](https://clades.nextstrain.org/dataset) HIV reference which states:

```
This data set uses the NCBI reference sequence NC_001802 based on the HXB2 genome K03455. The primary reason for choosing it is to ensure amino acid substitutions in conserved proteins such as Pol are numbered consistently. Note that this sequence as a few problems, including a premature stop-codon in nef.
```

-   `NC_001802.1.fasta`: Human immunodeficiency virus 1 genome fasta file downloaded from [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/NC_001802.1/)
-   `NC_001802.1.gff`: Human immunodeficiency virus 1 genome GFF3 annotation file downloaded from [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/NC_001802.1/)

#### `codfreq`

This reference was generated using the [HIV JSON profile](https://github.com/hivdb/codfreq/blob/main/profiles/HIV1.json) from codfreq software.

- `codfreq.fasta`: Was generated from the `"refSequence"` key of the .json file.
- `codonfreq.gff`: Was manually generated using the information from `"fragmentName"` and `"refRanges"` from `"fragmentConfig"`.

This is the default reference used in the nf-core/viralrecon HIV resistance detection protocol for the resulting codon frequencies and codon coverages to be directly comparable to those produced by [**HIVdb**](https://hivdb.stanford.edu/hivdb/by-reads/), ensuring accurate interpretation of resistance data.

### `illumina/sispa/`

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

### `illumina/amplicon`

| file                | num_seqs | sum_len   | min_len | avg_len | max_len | file_size | Sequencer   | LibrarySource |
| ------------------- | -------- | --------- | ------- | ------- | ------- | --------- | ----------- | ------------- |
| sample1_R1.fastq.gz | 27,721   | 8,285,732 | 35      | 168     | 301     | 4M        | PE Illumina | Metagenomics  |
| sample1_R2.fastq.gz | 27,721   | 8,285,900 | 35      | 168     | 301     | 4M        | PE Illumina | Metagenomics  |
| sample2_R1.fastq.gz | 21,481   | 6,416,734 | 35      | 168     | 301     | 3M        | PE Illumina | Metagenomics  |
| sample2_R2.fastq.gz | 21,481   | 6,416,265 | 35      | 168     | 301     | 3M        | PE Illumina | Metagenomics  |

> All FastQ files were sub-sampled to 0.02% of the original reads.

### `illumina/hiv/`

This dasatet was chosen because it is the example data for [HIVdb Drug Resistance Database](https://hivdb.stanford.edu/hivdb/by-reads/):

- DRR030302: Amplicon Whole Genome sequencing
- SRR4071760: Amplification of protease-RT genes
- SRR6937100: Amplification integrase gene

| file                  | num_seqs | sum_len   | min_len | avg_len | max_len | file_size | Sequencer   | LibrarySource |
| --------------------- | -------- | --------- | ------- | ------- | ------- | --------- | ----------- | ------------- |
| DRR030302_1.fastq.gz  | 10,512   | 2,545,102 | 40      | 242     | 251     | 1.3M      | PE Illumina | Viral RNA     |
| DRR030302_1.fastq.gz  | 10,512   | 2,545,205 | 40      | 242     | 251     | 1.8M      | PE Illumina | Viral RNA     |
| SRR4071760_1.fastq.gz | 10,582   | 2,524,863 | 45      | 238     | 251     | 1M        | PE Illumina | Synthetic     |
| SRR4071760_2.fastq.gz | 10,582   | 2,525,284 | 45      | 238     | 251     | 1.4M      | PE Illumina | Synthetic     |
| SRR6937100_1.fastq.gz | 10,484   | 1,295,077 | 35      | 123     | 151     | 556K      | PE Illumina | Genomic       |
| SRR6937100_2.fastq.gz | 10,484   | 1,289,631 | 33      | 123     | 151     | 612K      | PE Illumina | Genomic       |

> Original FastQ files were sub-sampled as explained in [Sampling procedure](#sampling-procedure)

### `illumina/fragmented/`

TBD

## Sampling procedure

### SARS-CoV-2

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


### HIV

The data was downsampled after Human Genome reads removal using different proportions:

- DRR030302: 0.025
- SRR4071760: 0.11
- SRR6937100: 0.55

We used the following commands:

```bash
seqtk sample -s100 <reads> <proportion>
```

## Expected output

TBD.
