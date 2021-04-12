# test-datasets: `rnaseq`

This branch contains test data to be used for automated testing with the [nf-core/rnaseq](https://github.com/nf-core/rnaseq) pipeline.

## Content of this repository

`reference/`: Sub-sampled genome reference files (iGenomes **S. cerevisiae** R64-1-1 Ensembl release)   

`testdata/*.fastq.gz`: Historical single-end test data for pipeline sub-sampled to ~2000 reads
`testdata/GSE110004/*.fastq.gz`: Paired-end test data for pipeline sub-sampled to 50000 reads

`samplesheet/samplesheet.csv`: Experiment design file for minimal test dataset  
`samplesheet/samplesheet_full.csv`: Experiment design file for full test dataset  

## Minimal test dataset origin

*S. cerevisiae* 101bp paired-end strand-specific RNA-seq dataset was obtained from:

> Andrew C K Wu, Harshil Patel, Minghao Chia, Fabien Moretto, David Frith, Ambrosius P Snijders, Folkert J van Werven. Repression of Divergent Noncoding Transcription by a Sequence-Specific Transcription Factor. Mol Cell. 2018 Dec 20;72(6):942-954.e7. doi: 10.1016/j.molcel.2018.10.018. [Pubmed](https://pubmed.ncbi.nlm.nih.gov/30576656/) [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110004)

### Sampling information

| run_accession | experiment_alias | read_count | sample_title                                                              |
|---------------|------------------|------------|---------------------------------------------------------------------------|
| SRR6357070    | GSM2879618       | 47629288   | Wild-type total RNA-Seq biological replicate 1                            |
| SRR6357071    | GSM2879619       | 68628914   | Wild-type total RNA-Seq biological replicate 2                            |
| SRR6357072    | GSM2879620       | 54771596   | Wild-type total RNA-Seq biological replicate 3                            |
| SRR6357073    | GSM2879621       | 56006930   | Rap1-AID degron no induction total RNA-Seq biological replicate 1         |
| SRR6357074    | GSM2879622       | 56259979   | Rap1-AID degron no induction total RNA-Seq biological replicate 2         |
| SRR6357075    | GSM2879623       | 51876040   | Rap1-AID degron no induction total RNA-Seq biological replicate 3         |
| SRR6357076    | GSM2879624       | 54935434   | Rap1-AID degron induction 30 minutes total RNA-Seq biological replicate 1 |
| SRR6357077    | GSM2879625       | 57770345   | Rap1-AID degron induction 30 minutes total RNA-Seq biological replicate 2 |
| SRR6357078    | GSM2879626       | 47537967   | Rap1-AID degron induction 30 minutes total RNA-Seq biological replicate 3 |
| SRR6357079    | GSM2879627       | 56870378   | Rap1-AID degron induction 2 hours total RNA-Seq biological replicate 1    |
| SRR6357080    | GSM2879628       | 59113530   | Rap1-AID degron induction 2 hours total RNA-Seq biological replicate 2    |
| SRR6357081    | GSM2879629       | 48202638   | Rap1-AID degron induction 2 hours total RNA-Seq biological replicate 3    |

### Sampling procedure

1. If we have a file called "chrI.fa" containing a single chromosome from _S. cerevisiae_ just edit the fasta entry header to include the taxonomy info as suggested in the Kraken2 manual (see [docs](https://github.com/DerrickWood/kraken2/wiki/Manual#custom-databases)) e.g. rename the entry header from `>I` to `>I|kraken:taxid|4932`.

    > NB: May not have to do this step but I just did it anyway.

2. Build Kraken2 database for custom genome

    ```console
    DBNAME='yeast_chrI'
    kraken2-build --download-taxonomy --db $DBNAME
    kraken2-build --add-to-library chrI.fa --db $DBNAME
    kraken2-build --build --db $DBNAME
    ```

3. (OPTIONAL) Download publicly available fastq files with [nf-core/rnaseq](https://github.com/nf-core/rnaseq) pipeline (see [docs](https://nf-co.re/rnaseq/3.0/usage#direct-download-of-public-repository-data)). This also auto-generates a samplesheet that can be easily re-formatted to work as input with nf-core/viralrecon in the following next step:

    ```console
    nextflow run nf-core/rnaseq \
        --public_data_ids ids.txt \
        -profile singularity
    ```

4. Only run the Kraken2 process from the [nf-core/viralrecon](https://github.com/nf-core/viralrecon) pipeline to get filtered fastq files:

    ```console
    nextflow run nf-core/viralrecon \
        --input samplesheet.csv \
        --kraken2_db yeast_chrI/ \
        --fasta chrI.fa \
        --platform illumina \
        --protocol metagenomic \
        --skip_fastqc \
        --skip_fastp \
        --skip_multiqc \
        --skip_assembly \
        --skip_variants \
        -profile singularity \
        -c custom.config \
    ```

    The contents of `custom.config` are defined below and are used to tell the pipeline to also publish the "fastq.gz" files because this isn't done by default:

    ```nextflow
    params {
        modules {
            'illumina_kraken2_run' {
                publish_files = ['txt':'', 'fastq.gz':'']
            }
        }
    }
    ```

5. The example command below was used to sub-sample the raw paired-end FastQ files to 50,000 reads (see [seqtk](https://github.com/lh3/seqtk)):

    ```console
    seqtk sample -s100 SRR6357070.classified_1.fastq.gz 50000 | gzip > SRR6357070_1.fastq.gz
    seqtk sample -s100 SRR6357070.classified_2.fastq.gz 50000 | gzip > SRR6357070_2.fastq.gz
    ```

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

## Create gff from gtf

In case the GTF gene annotation file gets updated, then GFF would also need to get updated. One can use [gffread](https://bioconda.github.io/recipes/gffread/README.html) to perform the conversion:

```console
gffread -F --keep-exon-attrs genes.gtf > genes.gff
```

Explanation of flags:

- `-F` preserves attributes for genes and transcripts, but doesn't preserve for exon features
- `--keep-exon-attrs` is needed as [featureCounts](http://subread.sourceforge.net/) in the [nf-core/rnaseq](https://github.com/nf-core/rnaseq/) pipeline uses the gene type/biotype (e.g. `protein_coding`, `lncRNA`) of the exons to count number of reads per biotype

## Create the gzipped references

In case the reference genomes or gene annotations get updated, the gzipped references would need to get updated, too. To make the gzipped references, run the following snippet in the `reference` folder:

```console
for F in $(ls -1 | grep -vE '.gz$'); do echo $F ; gzip -c $F > $F.gz ; done
```

This looks for files that don't end in `.gz` and compresses them.
