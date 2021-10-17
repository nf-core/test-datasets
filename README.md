# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

Test data to be used for automated testing with the nf-core pipelines

## Introduction

This branch contains test data for the [nf-core/nascent](https://github.com/nf-core/nascent) pipeline.

## Minimal test dataset origin

### Sampling information

### Sampling procedure

1. Downloaded GRCh38_chr21.fa from [NCBI](https://www.ncbi.nlm.nih.gov/nucleotide/CM000683.2). Click on "Send to:" and select "file"
2. Added >GRCh38_chr21|kraken:taxon|9606
3. Built the kraken db

```sh
DBNAME='GRCh38_chr21'
kraken2-build --download-taxonomy --db $DBNAME
kraken2-build --add-to-library chrI.fa --db $DBNAME
kraken2-build --build --db $DBNAME
```

4. Downloaded the fastqs using this script

```txt
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR014/SRR014285/SRR014285.fastq.gz -o SRR014285_Other_Sequencing_of_human.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR014/SRR014283/SRR014283.fastq.gz -o SRR014283_Other_Sequencing_of_human.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR014/SRR014284/SRR014284.fastq.gz -o SRR014284_Other_Sequencing_of_human.fastq.gz
```

```sh
nextflow run nf-core/viralrecon \
    --input results/samplesheet/samplesheet.csv \
    --kraken2_db GRCh38_chr21/ \
    --fasta GRCh38_chr21.fa \
    --platform illumina \
    --protocol metagenomic \
    --skip_fastqc \
    --skip_fastp \
    --skip_multiqc \
    --skip_assembly \
    --skip_variants \
    -profile singularity \
    -c custom.config -c ~/.nexflow.config
```

5. The example command below was used to sub-sample the raw paired-end FastQ files to 50,000 reads (see [rasusa](https://github.com/mbhall88/rasusa)):

```sh
rasusa -i results/kraken2/SRX882904_T2.classified.fastq.gz -c 50 -g 50000 -s 1 -o SRX882904_T2.fastq.gz
```

Or quickly for all the samples:

```sh
for i in results/kraken2/*.classified.fastq.gz; do
	echo "processing $i"
    filename=`basename "$i"`
    rasusa -i $i -c 50 -g 50000 -s 1 -o rasusa/$filename
done
```
