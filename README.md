# test-datasets: `rnaseq`
Test data to be used for automated testing with the nf-core pipelines

This branch contains test data for the [nf-core/rnaseq](https://github.com/nf-core/rnaseq) pipeline.

## Create gff from gtf

In case the GTF gene annotation file gets updated, then GFF would also need to get updated. One can use [gffread](https://bioconda.github.io/recipes/gffread/README.html) to perform the conversion:

```bash
gffread -F --keep-exon-attrs genes.gtf > genes.gff
```

Explanation of flags:

- `-F` preserves attributes for genes and transcripts, but doesn't preserve for exon features
- `--keep-exon-attrs` is needed as [featureCounts](http://subread.sourceforge.net/) in the [nf-core/rnaseq](https://github.com/nf-core/rnaseq/) pipeline uses the gene type/biotype (e.g. `protein_coding`, `lncRNA`) of the exons to count number of reads per biotype
