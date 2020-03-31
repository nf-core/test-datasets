# test-datasets: `covid19`

This branch contains test data to be used for automated testing with the [nf-core/covid19](https://github.com/nf-core/covid19) pipeline.

## Content of this repository

### `samplesheet_SISPA.csv`

Sample information sheet required to test the pipeline containing sample information and links to FastQ files stored in this repository. This sample sheet corresponds to Illumina SISPA data.

### `reference/`

`hg19_chr21.fa`: Entire chromosome 21 DNA sequence from the UCSC hg19 human assembly.
`NC_045512.2.fasta`: Reference SARS-Cov2 fasta file.
`NC_045512.2.gff`: Reference SARS-Cov2 gff file.

### `illumina/`

| file                    | num_seqs | sum_len    | min_len | avg_len | max_len | file_size | Sequencer   | LibrarySource      |
|-------------------------|----------|------------|---------|---------|---------|-----------|-------------|--------------------|
| SRR10903401_1.fastq.gz  |    9,580 |  1,442,283 |     143 |   150.6 |     151 |      665K | PE Illumina | Metatranscriptomic |
| SRR10903401_2.fastq.gz  |    9,580 |  1,443,343 |     145 |   150.7 |     151 |      751K | PE Illumina | Metatranscriptomic |
| SRR10903402_1.fastq.gz  |   13,588 |  2,045,955 |     144 |   150.6 |     151 |      991K | PE Illumina | Metatranscriptomic |
| SRR10903402_2.fastq.gz  |   13,588 |  2,046,944 |     141 |   150.6 |     151 |      1.2M | PE Illumina | Metatranscriptomic |
| SRR11092056_1.fastq.gz  |  104,641 | 14,823,897 |      35 |   141.7 |     151 |      9.7M | PE Illumina | Metagenomics       |
| SRR11092056_2.fastq.gz  |  104,641 | 14,823,709 |      35 |   141.7 |     151 |      11M  | PE Illumina | Metagenomics       |
| SRR11177792_1.fastq.gz  |  104,780 | 28,110,508 |      35 |   268.3 |     301 |      15M  | PE Illumina | Genomic            |
| SRR11177792_2.fastq.gz  |  104,780 | 28,110,045 |      35 |   268.3 |     301 |      16M  | PE Illumina | Genomic            |
| SRR11241255.fastq.gz    |    2,442 |    438,067 |      20 |   179.4 |     185 |      86K  | SE Illumina | Viral RNA          |
| SRR11140744_R1.fastq.gz |    10092 |    2284737 |     100 |   175.5 |     251 |      747K | PE Illumina | Metagenomics       |
| SRR11140744_R2.fastq.gz |    10092 |    2260970 |     100 |   175.5 |     251 |      783K | PE Illumina | Metagenomics       |
| SRR11140746_R1.fastq.gz |     7196 |    1609884 |     100 |   175.5 |     251 |      554K | PE Illumina | Metagenomics       |
| SRR11140746_R2.fastq.gz |     7196 |    1594703 |     100 |   175.5 |     251 |      580K | PE Illumina | Metagenomics       |
| SRR11140748_R1.fastq.gz |     8447 |    1918541 |     100 |   175.5 |     251 |      650K | PE Illumina | Metagenomics       |
| SRR11140748_R2.fastq.gz |     8447 |    1903781 |     100 |   175.5 |     251 |      683K | PE Illumina | Metagenomics       |
| SRR11140750_R1.fastq.gz |      369 |      81898 |     100 |   175.5 |     251 |       40K | PE Illumina | Metagenomics       |
| SRR11140750_R2.fastq.gz |      369 |      80344 |     102 |   176.5 |     251 |       41K | PE Illumina | Metagenomics       |


> All FastQ files were sub-sampled to 0.02% of the original reads.

### `nanopore/`

| file                    | num_seqs | sum_len    | min_len | avg_len | max_len | file_size | Sequencer   | LibrarySource      |
|-------------------------|----------|------------|---------|---------|---------|-----------|-------------|--------------------|
| SRR10948474.fastq.gz    |   10,141 |  5,655,211 |     125 |   557.7 |   5,892 |      5.7M |    Nanopore | Genomic            |
| SRR10948550.fastq.gz    |    8,539 |  2,953,453 |     104 |   345.9 |   1,577 |      2.9M |    Nanopore | Genomic            |

> All FastQ files were sub-sampled to 0.02% of the original reads.

## Sampling procedure

Prepare a file `list.txt` with the following SRA accession numbers:

```
SRR10903401
SRR10948474
SRR10903402
SRR10948550
SRR11092056
SRR11177792
SRR11241255
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
