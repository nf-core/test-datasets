# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the nf-core pipelines

This branch contains test data to be used for automated testing with the [nf-core/fastquorum](https://github.com/nf-core/fastquorum) pipeline.

## Content of this repository

`references/`: genome reference and auxiliary files for Homo sapiens assembly hg38 for chromosome 17.

`testdata/duplex-seq`: raw FASTQs for Illumina paired-end duplex-sequencing experiments


The full contents are shown below:

```console
.
├── CITATION.cff
├── LICENSE
├── README.md
├── docs
│   ├── ADD_NEW_DATA.md
│   ├── USE_EXISTING_DATA.md
│   └── images
│       ├── test-datasets_logo.png
│       └── test-datasets_logo.svg
├── references
│   ├── chr17.dict
│   ├── chr17.fa
│   ├── chr17.fa.amb
│   ├── chr17.fa.ann
│   ├── chr17.fa.bwt
│   ├── chr17.fa.fai
│   ├── chr17.fa.pac
│   ├── chr17.fa.sa
│   └── samplesheet.csv
└── testdata
    └── duplex-seq
        ├── SRR6109255_1.fastq.gz
        ├── SRR6109255_2.fastq.gz
        ├── SRR6109273_1.fastq.gz
        └── SRR6109273_2.fastq.gz
```

### Sample Information

| Run Accession | Experiment Accession | Experiment Title                                        | Citation |
|---------------|----------------------|---------------------------------------------------------|----------|
| SRR6109255    | SRX3224128           | Illumina MiSeq sequencing: CRISPR-DS Sequencing of TP53 | [1]      |
| SRR6109273    | SRX3224110           | Illumina MiSeq sequencing: CRISPR-DS Sequencing of TP53 | [1]      |

Citations:

1. Nachmanson D, Lian S, Schmidt EK, Hipp MJ, Baker KT, Zhang Y, Tretiakova M, Loubet-Senear K, Kohrn BF, Salk JJ, Kennedy SR, Risques RA. Targeted genome fragmentation with CRISPR/Cas9 enables fast and efficient enrichment of small genomic regions and ultra-accurate sequencing with low DNA input (CRISPR-DS). Genome Res. 2018 Oct;28(10):1589-1599. doi: 10.1101/gr.235291.118. Epub 2018 Sep 19. PMID: 30232196; PMCID: PMC6169890.

For `SRR6109255`, we sub-sampled the reads to fit under the GitHub 100MB limit:

```console
seqtk sample -s 42 SRR6109255_1.fastq.gz 0.65 | gzip -c > full/SRR6109255_1.fastq.gz
seqtk sample -s 42 SRR6109255_2.fastq.gz 0.65 | gzip -c > full/SRR6109255_2.fastq.gz
```

We also sub-sampled `SRR6109255` to create a tiny dataset for rapid testing:

```console
seqtk sample -s 42 SRR6109255_1.fastq.gz 0.01 | gzip -c > tiny/SRR6109255_1.fastq.gz
seqtk sample -s 42 SRR6109255_2.fastq.gz 0.01 | gzip -c > tiny/SRR6109255_2.fastq.gz
```


|    tool |  version |
|---------|----------|
| `seqtk` | 1.4-r122 |

### Reference Information

* `chr17` downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr17.fa.gz.

We used the following commands to prepare auxiliary reference information:

```console
samtools faidx chr17.fa
samtools dict -u https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr17.fa.gz -a hg38 -s "Homo sapiens" chr17.fa > chr17.dict
bwa index chr17.fa
```
|       tool |      version |
|------------|--------------|
| `samtools` |         1.17 |
|      `bwa` | 0.7.17-r1188 |

