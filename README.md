# test-datasets: `viralrecon`

This branch contains test data to be used for automated testing with the [nf-core/viralrecon](https://github.com/nf-core/viralrecon) pipeline.

## Content of this repository

### `samplesheet.csv`

Sample information sheet required to test the pipeline containing sample information and links to FastQ files stored in this repository.

### `reference/`

`hg19_chr21.fa`: Entire chromosome 21 DNA sequence from the UCSC hg19 human assembly.

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

> All FastQ files were sub-sampled to 0.02% of the original reads.

## Sampling procedure

Prepare a file `list.txt` with the following SRA accession numbers:

```
SRR10903401
SRR10903402
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
