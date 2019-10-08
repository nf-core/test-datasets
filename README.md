# test-datasets: `rnaseq`
Test data to be used for automated testing with the nf-core pipelines

This branch contains test data for the [nf-core/rnaseq](https://github.com/nf-core/rnaseq) pipeline.

## Reference indices

### HISAT2 indices

HiSAT2 indices were made from the reduced chromosome I data from the nf-core/rnaseq test dataset, by using `--saveReference`, and then using `tar` on them:

```
cd results/reference_genome
tar -zcvf hisat2.tar.gz *.hisat2_*
```

### STAR index

STAR indices were made from the reduced chromosome I data from the nf-core/rnaseq test dataset, by using `--saveReference`, and then using `tar` on them:

```
cd results/reference_genome
tar -zcvf star.tar.gz star
```

### Salmon index

Salmon indices were made from the reduced chromosome I data from the nf-core/rnaseq test dataset, by using `--saveReference`, and then using `tar` on them:

```
cd results/reference_genome
tar -zcvf salmon_index.tar.gz salmon_index
```
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
