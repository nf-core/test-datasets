# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the nf-core pipelines

## Introduction

nf-core is a collection of high quality Nextflow pipelines.

## Documentation
nf-core/test-datasets comes with documentation in the `docs/` directory:

01. [Add a new  test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/ADD_NEW_DATA.md)
02. [Use an existing test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/USE_EXISTING_DATA.md)

## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).

## Datasets for nf-core/magmap

### Data added to implement genome preference

This is data created to deal with [issue #187](https://github.com/nf-core/magmap/issues/187).
Used in the test config `test_sourmash_mix_dupl_species`.

Seven archaeal species selected from the Sulfolobales order, each having two genomes available at GTDB (R10-RS226):

d__Archaea;p__Thermoproteota;c__Thermoprotei_A;o__Sulfolobales;f__Acidilobaceae;g__Aeropyrum;s__Aeropyrum pernix
d__Archaea;p__Thermoproteota;c__Thermoprotei_A;o__Sulfolobales;f__Ignicoccaceae;g__Ignicoccus;s__Ignicoccus hospitalis
d__Archaea;p__Thermoproteota;c__Thermoprotei_A;o__Sulfolobales;f__Ignisphaeraceae;g__Ignisphaera;s__Ignisphaera cupida
d__Archaea;p__Thermoproteota;c__Thermoprotei_A;o__Sulfolobales;f__Pyrodictiaceae;g__Pyrodictium;s__Pyrodictium delaneyi
d__Archaea;p__Thermoproteota;c__Thermoprotei_A;o__Sulfolobales;f__Sulfolobaceae;g__Metallosphaera;s__Metallosphaera hakonensis
d__Archaea;p__Thermoproteota;c__Thermoprotei_A;o__Sulfolobales;f__Sulfolobaceae;g__Metallosphaera;s__Metallosphaera javensis
d__Archaea;p__Thermoproteota;c__Thermoprotei_A;o__Sulfolobales;f__Sulfolobaceae;g__Saccharolobus;s__Saccharolobus caldissimus

Contigs for all 14 genomes (both for each species) were downloaded from NCBI and were used to:

* Generate 1000 read pairs per genome

* The 14 read pair files were concatenated 5, 5 and 4 into sample files:

testdata/archaeal_duplicates00_1.fastq.gz
testdata/archaeal_duplicates00_2.fastq.gz
testdata/archaeal_duplicates01_1.fastq.gz
testdata/archaeal_duplicates01_2.fastq.gz
testdata/archaeal_duplicates02_1.fastq.gz
testdata/archaeal_duplicates02_2.fastq.gz

* The three read pair files where included in a new samplesheet:

samplesheets/archaeal_duplicate_genomes_per_species.csv

* The species-representative genomes were used to generate a Sourmash index (archaeal_duplicates.index.sbt.zip)

testdata/archaeal_duplicates.index.sbt.zip

* A genomeinfo csv plus CheckM, CheckM2 and GTDB-Tk files were created for the non species-representative from the official GTDB metadata

testdata/archaeal_duplicates.genomes.csv
testdata/archaeal_duplicates.checkm2.tsv
testdata/archaeal_duplicates.checkm.tsv

