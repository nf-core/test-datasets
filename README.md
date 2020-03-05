# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

Test data to be used for automated testing with the nf-core pipelines

## Test data and references for Sarek

This branch contains test data and references for the [nf-core/sarek](https://github.com/nf-core/sarek) pipeline.

### Test data

Small files sub-sampled and restricted for testing

* `target.bed`

For targeted sequencing.

* `dummy/`, `manta/` and `tiny/`

> Used to test `Preprocessing` only or in conjunction with following steps.

Contain `normal/` and `tumor/` fastq paired files for sarek.

* `recalbam/`

> Used to test `Variant Calling`.

Contain `bam` and indexes made from `manta/` and `dummy/` fastq files using `Sarek`.

* `vcf/`

> Used to test `Annotation`.

Contain VCF files made with `Strelka` using `Sarek`.

* `tsv/`

Contain TSV files used to launch `Sarek` with different settings.

* `*-https.tsv` files use `HTTPS` paths to files, and the other files will use local relative path.
* `*-normal*` files contain only `normal` sample, and the other files will use `normal` and `tumor` sample.
* `*-recal*` files contain bam files and are to be used for `Variant Calling`.
* `*-multiple*` files contain multiple samples.



* `umi-dna`

Used to test UMI functionality. Contains two different folders:

- Illumina TSO

This small dataset includes only reads mapped to chromosome 20. A smaller reference could also be used.
Both reads include UMIs with structure 7M1S+T (i.e. a 7bp UMI followed by a constant base and the template)

- QIAseq reads

From SRA repository, simply selected the first 300,000 reads from SRR7545951: this sample has been used to develop smCounter2.
In this case only read 2 includes a UMI, with structure 12M11S+T.

### Reference

Small reference files for testing, based on `GRCh37`

* `1000G_phase1.indels.b37.small.vcf.gz`
  Indels and index

* `1000G_phase3_20130502_SNP_maf0.3.small.loci`, `1000G_phase3_20130502_SNP_maf0.3.small.loci.gc`
  Files for `ASCAT`

* `dbsnp_138.b37.small.vcf.gz`
  DBSNP and index

* `human_g1k_v37_decoy.small.fasta`
  Reference genome and indexes

* `Mills_and_1000G_gold_standard.indels.b37.small.vcf.gz`
  Indels and index

* `small.intervals`
  Intervals for Preprocessing and Variant Calling

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).




