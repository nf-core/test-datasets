# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the nf-core pipelines

## Introduction

nf-core is a collection of high quality Nextflow pipelines. This repository contains various files for CI and unit testing of nf-core pipelines and infrastructure.

The principle for nf-core test data is as small as possible, as large as necessary. Always ask for guidance on the [nf-core slack](https://nf-co.re/join) before adding new test data.

## Documentation

nf-core/test-datasets comes with documentation in the `docs/` directory:

01. [Add a new  test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/ADD_NEW_DATA.md)
02. [Use an existing test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/USE_EXISTING_DATA.md)

## Downloading test data

Due the large number of large files in this repository for each pipeline, we highly recommend cloning only the branches you would use.

```bash
git clone <url> --single-branch --branch <pipeline/modules/branch_name>
```

To subsequently clone other branches[^1]

```bash
git remote set-branches --add origin [remote-branch]
git fetch
```

## nf-core/taxprofiler specific information

### fastq


The main test data used for nf-core/taxprofiler is from [Maixner et al. (2021) _Curr. Bio._](https://doi.org/10.1016/j.cub.2021.09.031), with ENA project accession ID: PRJEB44507. The following selected libraries were all sequenced on an Illumina MiSeq, and were selected due to their small size (~1million reads, <100MB) and known mixture of (gut) bacteria, (ancient human) eukaryotes, and (yeast) fungi (according to the results of the paper).

- ERX5474937
- ERX5474932
- ERX5474930
- ERX5474936

Data was downloaded with nf-core/fetchNGS 1.5 (with Nextflow 21.10.06):

```bash
nextflow run nf-core/fetchngs --input maixner2021_acc_codes.txt --input_type sra
```
FASTQ files are stored under `data/fastq/`

Test data for long reads are the same as [nf-core/mag](https://github.com/nf-core/mag) is using and can be found [here](https://github.com/nf-core/test-datasets/blob/mag/test_data/minigut_reads.fastq.gz)
### fasta

One of the files was converted to FASTA file with seqtk 1.3-r106

```bash
seqtk seq -a  ERX5474930_ERR5766174_1.fastq.gz > ERX5474930_ERR5766174_1.fa.gz 
```

FASTA files are stored under `data/fasta/`

### databases

Two abundant species found in the above dataset are _Penicillium roqueforti_
and _Saccharomyces cerevisiae_. The genomes of the two P. roqueforti were 
downloaded with:

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/533/775/GCF_015533775.1_ASM1553377v1/GCF_015533775.1_ASM1553377v1_genomic.fna.gz # P. roqueforti
```

In addition we include the human mitochondrial genome.

```bash
curl 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=251831106&extrafeat=null&conwithfeat=on&hide-cdd=on'| gzip > NC_012920.1.fa.gz # H. sapiens mito
```

#### MALT

```bash
wget https://software-ab.informatik.uni-tuebingen.de/download/megan6/megan-nucl-Feb2022.db.zip
unzip megan-nucl-Feb2022.db
malt-build -i *.gz -s DNA -d taxprofiler-testdb -t 8 -st 4 -a2t megan-nucl-Feb2022.db
```
#### Kraken2

```bash
kraken2-build --download-taxonomy --db taxprofiler-testdb
gunzip ../raw/*.gz
kraken2-build --add-to-library ../raw/GCF_000146045.2_R64_genomic.fna --db taxprofiler-testdb/
kraken2-build --add-to-library ../raw/NC_012920.1.fa --db taxprofiler-testdb/
kraken2-build --build --db taxprofiler-testdb/
kraken2-build --clean --db taxprofiler-testdb/
```
#### Centrifuge
We use the same test database for [centrifuge](https://github.com/nf-core/test-datasets/blob/mag/test_data/minigut_cf.tar.gz) as [nf-core/mag](https://github.com/nf-core/mag) is using in their test data.

## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).

[^1]: From [stackoverflow](https://stackoverflow.com/a/60846265/11502856)
