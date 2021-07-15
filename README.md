# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the nf-core pipelines

## nf-core/circrna
This branch contains test-data for `nf-core/circrna`.

### Contents of branch:
* `fastq/` 6 FASTQ read pairs.
* `reference/` Reference annotation files
* `phenotype.csv` metadata file for DESeq2.
* `samples.csv` input test-dataset csv file

### Test-dataset generation strategy:
Sequencing reads retrieved using `nf-core/fetchngs`, input file contained the accession `PRJNA669975`.

Reads were mapped to the `ce10` genome using hisat2. Bam files were subset according to a specific region:

```console
for bam in *.bam; do samtools view -b $bam chrI:1000000-4000000 > /data/bdigby/igenomes_test/bams/${bam%.bam}.bam; done
```

And subsequently converted to FASTQ files using `Picard`:

```console
for bam in *.bam; do picard -Xmx2g SamToFastq I=$bam F=${bam%.bam}_1.fastq.gz F2=${bam%.bam}_2.fastq.gz VALIDATION_STRINGENCY=LENIENT ; done
```

Reference FASTA and GTF file were subset to contain only `chrI`.

- Barry
