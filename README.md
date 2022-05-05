# test-datasets: `rnavar`

This branch contains test data to be used for automated testing with the [nf-core/rnavar](https://github.com/nf-core/rnavar) pipeline.

## Content of this repository

`samplesheet/samplesheet.csv`: Experiment design file for minimal test dataset
`samplesheet/samplesheet_full.csv`: Experiment design file for full test dataset

## Minimal test dataset origin

The raw data was retrieved from 'Genome in a Bottle' sample GM12878 (SRA accession [SRX2900878](https://www.ncbi.nlm.nih.gov/sra/?term=SRX2900878), sequenced on NextSeq 500 with 151bpx2 library).

### Sampling procedure

1. The data was downloaded using the SRA Toolkit with:

    ```bash
    prefetch -v <Acc number>
    sam-dump <Acc number> | samtools view -bS - > <Acc number>.bam
    ```

2. Manually verified a region covered by reads in chromsome 22 using IGV genome veiwer.

    ```bash
    chr22   16570000        16610000
    ```

3. Extracted all the reads and their pair-mates that overlap the SNP sites in the above region were extracted using [VariantBAM](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4920121/) and converted to `fastq.gz` using [qbic-pipelines/bamtofastq](https://github.com/qbic-pipelines/bamtofastq).

## Full test dataset origin

*H. sapiens* paired-end strand-specific RNA-seq dataset was obtained from:

> ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the human genome. Nature 2012 Sep 6;489(7414):57-74. [Pubmed](https://pubmed.ncbi.nlm.nih.gov/22955616/)

The GM12878 and K562 ENCODE data was also used to benchmark RNA-seq quantification pipelines in the paper below:

> Mingxiang Teng, Michael I. Love, Carrie A. Davis, Sarah Djebali, Alexander Dobin, Brenton R. Graveley, Sheng Li, Christopher E. Mason, Sara Olson, Dmitri Pervouchine, Cricket A. Sloan, Xintao Wei, Lijun Zhan, and Rafael A. Irizarry. A benchmark for RNA-seq quantification pipelines. Genome Biol. 2016; 17: 74. Published online 2016 Apr 23. doi: 10.1186/s13059-016-0940-1. [Pubmed](https://pubmed.ncbi.nlm.nih.gov/27107712/)


| study_alias | run_accession | experiment_alias | encode_library_id | sample_description | instrument_model | library_layout | read_count | sex | fastq_ftp | fastq_md5 |
|-------------|---------------|------------------|-------------------|--------------------|------------------|----------------|------------|-----|-----------|-----------|
  | [GSE78551](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78551) | SRR3192657 | GSM2072350 | ENCLB038ZZZ | Homo sapiens GM12878 immortalized cell line | Illumina HiSeq 2000 | PAIRED | 93555584 | female | [fastq_1](ftp.sra.ebi.ac.uk/vol1/fastq/SRR319/007/SRR3192657/SRR3192657_1.fastq.gz) [fastq_2](ftp.sra.ebi.ac.uk/vol1/fastq/SRR319/007/SRR3192657/SRR3192657_2.fastq.gz) | f3a3aee0e1f0f54dc9afd8f7c0442aba;6bff7e7d944736251cfbc36e35c3f431 |
| [GSE78551](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78551) | SRR3192658 | GSM2072351 | ENCLB037ZZZ | Homo sapiens GM12878 immortalized cell line | Illumina HiSeq 2000 | PAIRED | 97548052 | female | [fastq_1](ftp.sra.ebi.ac.uk/vol1/fastq/SRR319/008/SRR3192658/SRR3192658_1.fastq.gz) [fastq_2](ftp.sra.ebi.ac.uk/vol1/fastq/SRR319/008/SRR3192658/SRR3192658_2.fastq.gz) | f6fdb08100033d98bfcba0801a838bf9;b369f63c5d37e515b4e102fa8c8d75e7 |
| [GSE78557](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78557) | SRR3192408 | GSM2072362 | ENCLB055ZZZ | Homo sapiens K562 immortalized cell line | Illumina HiSeq 2000 | PAIRED | 92172367 | female | [fastq_1](ftp.sra.ebi.ac.uk/vol1/fastq/SRR319/008/SRR3192408/SRR3192408_1.fastq.gz) [fastq_2](ftp.sra.ebi.ac.uk/vol1/fastq/SRR319/008/SRR3192408/SRR3192408_2.fastq.gz) | 53815dcaeeb331459ab72bffe0a9432f;e73d0e7b764d96f08cf2caf4a7e880ff |
| [GSE78557](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78557) | SRR3192409 | GSM2072363 | ENCLB056ZZZ | Homo sapiens K562 immortalized cell line | Illumina HiSeq 2000 | PAIRED | 113327735 | female | [fastq_1](ftp.sra.ebi.ac.uk/vol1/fastq/SRR319/009/SRR3192409/SRR3192409_1.fastq.gz) [fastq_2](ftp.sra.ebi.ac.uk/vol1/fastq/SRR319/009/SRR3192409/SRR3192409_2.fastq.gz) | 5904c8781f4fd6771a5e9a32696cd49b;b23e23639258c93944ff9a64b08b9f67 |
| [GSE90237](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90237) | SRR5048099 | GSM2400174 | ENCLB555AQN | Homo sapiens MCF-7 immortalized cell line | Illumina Genome Analyzer IIx | PAIRED | 128178110 | female | [fastq_1](ftp.sra.ebi.ac.uk/vol1/fastq/SRR504/009/SRR5048099/SRR5048099_1.fastq.gz) [fastq_2](ftp.sra.ebi.ac.uk/vol1/fastq/SRR504/009/SRR5048099/SRR5048099_2.fastq.gz) | c23adfcad78e9162a83e18fc76e7ebfd;fd0c3baabd67659aecf6c88feef30259 |
| [GSE90237](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90237) | SRR5048100 | GSM2400175 | ENCLB555AQO | Homo sapiens MCF-7 immortalized cell line | Illumina Genome Analyzer IIx | PAIRED | 131814222 | female | [fastq_1](ftp.sra.ebi.ac.uk/vol1/fastq/SRR504/000/SRR5048100/SRR5048100_1.fastq.gz) [fastq_2](ftp.sra.ebi.ac.uk/vol1/fastq/SRR504/000/SRR5048100/SRR5048100_2.fastq.gz) | f7e732c768e4080311a49e6048c4d515;5619f168e72c5ca27b1b805a91de4444 |
| [GSE90225](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90225) | SRR5048077 | GSM2400152 | ENCLB555AMA | Homo sapiens H1-hESC stem cell male embryo | Illumina Genome Analyzer IIx | PAIRED | 125395196 | male | [fastq_1](ftp.sra.ebi.ac.uk/vol1/fastq/SRR504/007/SRR5048077/SRR5048077_1.fastq.gz) [fastq_2](ftp.sra.ebi.ac.uk/vol1/fastq/SRR504/007/SRR5048077/SRR5048077_2.fastq.gz) | 6beb20b2cd99542433986b8fe844ef09;4f63ef9e16dc9f0f8be159b02d40f0c6 |
| [GSE90225](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90225) | SRR5048078 | GSM2400153 | ENCLB555AMB | Homo sapiens H1-hESC stem cell male embryo | Illumina Genome Analyzer IIx | PAIRED | 107101340 | male | [fastq_1](ftp.sra.ebi.ac.uk/vol1/fastq/SRR504/008/SRR5048078/SRR5048078_1.fastq.gz) [fastq_2](ftp.sra.ebi.ac.uk/vol1/fastq/SRR504/008/SRR5048078/SRR5048078_2.fastq.gz) | 9c60d407bae58019889b13acb1032116;fc5df7d28daf6df1b212aaac914f1324 |
