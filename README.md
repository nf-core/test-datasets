# test-datasets: `rnaseq`
Test data to be used for automated testing with the nf-core pipelines

This branch contains test data for the [nf-core/rnaseq](https://github.com/nf-core/rnaseq) pipeline.

## Create the gzipped references
In case the reference genomes or gene annotations get updated, the gzipped references would need to get updated, too. To make the gzipped references, run the following snippet in the `reference` folder:

```bash
for F in $(ls -1 | grep -vE '.gz$'); do echo $F ; gzip -c $F > $F.gz ; done
```

This looks for files that don't end in `.gz` and compresses them.
