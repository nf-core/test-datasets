# test-datasets: `viralrecon`

This branch contains test data to be used for automated testing with the [nf-core/viralrecon](https://github.com/nf-core/viralrecon) pipeline.

## Content of this repository

### `samplesheet/`

This directory contains the sample sheets used to test different test configurations of the `nf-core/viralrecon` pipeline.

#### `samplesheet_test_nanopore.csv`

Sample sheet for Nanopore test data.  
Includes SARS-CoV-2 sample names and MinION barcodes hosted in this repository.  
For testing purposes, some barcodes are included without associated data.  
Raw data can be found in [`nanopore/minion`](nanopore/minion) and have been subsetted to include a maximum of three `fast5` or `fastq` files per barcode.

#### `samplesheet_test_sra.csv`

Sample sheet for SISPA-based Illumina data.  
Contains SARS-CoV-2 SISPA probe enriched sample information, including one link to original full FastQ files and two datasets to be downloaded from SRA (one single-end and one paired-end).  
Used to test SRA-based inputs and mixed dataset handling.

#### `samplesheet_test_illumina_sispa.csv`

Sample sheet for small-scale SISPA test data.  
Contains SARS-CoV-2 SISPA probe enriched sample information and links to corresponding FastQ files.  

#### `samplesheet_test_illumina_amplicon.csv`

Sample sheet for small-scale amplicon test data.  
Includes SARS-CoV-2 amplicon primer enriched sample sample information and links to FastQ files hosted in this repository.  

#### `samplesheet_full_illumina_sispa.csv`

Sample sheet for full-scale SISPA test data.  
Contains SARS-CoV-2 SISPA probe enriched sample information and links to corresponding FastQ files.  

#### `samplesheet_full_illumina_amplicon.csv`

Sample sheet for full-scale amplicon test data.  
Includes SARS-CoV-2 amplicon primer enriched sample sample information and links to FastQ files hosted in this repository.  

#### `samplesheet_full_illumina_fragmented.csv`

Sample sheet for fragmented genome tests.  
Contains Crimean-Congo hemorrhagic fever virus sample information and links to FastQ files stored in this repository.  
Used to assess pipeline performance on non-contiguous viral genomes.

#### `v3.0/samplesheet_test_hiv.csv`

Sample sheet for HIV  test data.  
Contains HIV Illumina amplicon primer enriched sample information and links to FastQ files stored in this repository.  
Data originate from multiple SRA experiments and are treated as non-amplicon data for testing purposes.

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

Small custom blast database containing 8 enterovirus genomes with taxid mapping required to test the pipeline. The steps used to generate the database are:

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
makeblastdb -in genomes.fasta -dbtype nucl -parse_seqids -out minimal_ev_database -taxid_map taxid_map.txt
```

Manually download the files `taxdb.bt[id]` and `taxonomy4blast.sqlite3` from NCBI, through [this link](https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz/) and add to the created blast database folder

Compress using `tar`

```bash
tar -czvf minimal_ev_db.tar.gz minimal_ev_database/
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

#### `crimea_congo`

Crimea Congo is a fragmented genome with three fragments. S, M and L based on the fragment size.

-   `crimea_congo.fasta.gz`: Crimea Congo fasta genome containing S, M and L fragments: KY484036.1, KY484035.1, KY484034.1
-   `crimea_congo.gff.gz`: Crimea congo genome GFF3 annotation file containing annotation for S, M and L fragments: KY484036.1, KY484035.1, KY484034.1

#### `NC_001802.1`

This reference was chosen based on [Nextclade's](https://clades.nextstrain.org/dataset) HIV reference which states:

```
This data set uses the NCBI reference sequence NC_001802 based on the HXB2 genome K03455. The primary reason for choosing it is to ensure amino acid substitutions in conserved proteins such as Pol are numbered consistently. Note that this sequence has a few problems, including a premature stop-codon in nef.
```

-   `NC_001802.1.fasta`: Human immunodeficiency virus 1 genome fasta file downloaded from [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/NC_001802.1/)
-   `NC_001802.1.gff`: Human immunodeficiency virus 1 genome GFF3 annotation file downloaded from [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/NC_001802.1/)

#### `codfreq`

This reference was generated using the [HIV JSON profile](https://github.com/hivdb/codfreq/blob/main/profiles/HIV1.json) from [codfreq](https://github.com/hivdb/codfreq) software.

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

This dataset was chosen because it is the example data for [HIVdb Drug Resistance Database](https://hivdb.stanford.edu/hivdb/by-reads/):

- DRR030302: Amplicon Whole Genome sequencing
- SRR4071760: Amplification of protease-RT genes
- SRR6937100: Amplification of integrase genes

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

### `fastq/illumina_enterovirus/`

| file                   | num_seqs | sum_len   | min_len | avg_len | max_len | file_size | Sequencer   | LibrarySource |
| ---------------------- | -------- | --------- | ------- | ------- | ------- | --------- | ----------- | ------------- |
| SRR13266665_1.fastq.gz | 10,000   | 1,342,331 | 35      | 134.23  | 151     | 508K      | PE Illumina | Metagenomics  |
| SRR13266665_2.fastq.gz | 10,000   | 1,321,544 | 35      | 132.15  | 151     | 508K      | PE Illumina | Metagenomics  |

> All FastQ files were sub-sampled down to 10,000 reads.

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

### Enterovirus

The data was sub-sampled to 10,000 reads using `seqtk`

```bash
seqtk sample -s100 <read1> 10000 > <output>
```

## Expected output

TBD.
