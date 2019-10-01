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
