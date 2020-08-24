# test-datasets: `rnaseq`

This branch contains test data to be used for automated testing with the [nf-core/rnaseq](https://github.com/nf-core/rnaseq) pipeline.

## Content of this repository

`reference/`: Sub-sampled genome reference files (iGenomes **S. cerevisiae** R64-1-1 Ensembl release)   
`testdata/*.fastq.gz`: Historical single-end test data for pipeline sub-sampled to ~2000 reads
`testdata/GSE110004/*.fastq.gz`: Paired-end test data for pipeline sub-sampled to ~100000 reads
`samplesheet/samplesheet.csv`: Experiment design file for minimal test dataset  
`samplesheet/samplesheet_full.csv`: Experiment design file for full test dataset  

## Minimal test dataset origin

*S. cerevisiae* paired-end strand-specifc RNA-seq dataset was obtained from:

Andrew C K Wu, Harshil Patel, Minghao Chia, Fabien Moretto, David Frith, Ambrosius P Snijders, Folkert J van Werven. Repression of Divergent Noncoding Transcription by a Sequence-Specific Transcription Factor. Mol Cell. 2018 Dec 20;72(6):942-954.e7. doi: 10.1016/j.molcel.2018.10.018. [Pubmed](https://pubmed.ncbi.nlm.nih.gov/30576656/) [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110004)

### Sample information

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

The example command below was used to sub-sample the raw paired-end FastQ files to 100,000 reads (see [seqtk](https://github.com/lh3/seqtk)).

```bash
mkdir -p sample
seqtk sample -s100 SRR6357070_1.fastq.gz 100000 | gzip > ./sample/SRR6357070_1.fastq.gz
seqtk sample -s100 SRR6357070_2.fastq.gz 100000 | gzip > ./sample/SRR6357070_2.fastq.gz
```

## Full test dataset origin

TBD.

## Create gff from gtf

In case the GTF gene annotation file gets updated, then GFF would also need to get updated. One can use [gffread](https://bioconda.github.io/recipes/gffread/README.html) to perform the conversion:

```bash
gffread -F --keep-exon-attrs genes.gtf > genes.gff
```

Explanation of flags:

- `-F` preserves attributes for genes and transcripts, but doesn't preserve for exon features
- `--keep-exon-attrs` is needed as [featureCounts](http://subread.sourceforge.net/) in the [nf-core/rnaseq](https://github.com/nf-core/rnaseq/) pipeline uses the gene type/biotype (e.g. `protein_coding`, `lncRNA`) of the exons to count number of reads per biotype

## Create the gzipped references

In case the reference genomes or gene annotations get updated, the gzipped references would need to get updated, too. To make the gzipped references, run the following snippet in the `reference` folder:

```bash
for F in $(ls -1 | grep -vE '.gz$'); do echo $F ; gzip -c $F > $F.gz ; done
```

This looks for files that don't end in `.gz` and compresses them.
