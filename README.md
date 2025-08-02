# Test data for methylong pipeline

This branch contains test data for testing the pipeline `methylong` .

## Minimal Test data

The minimal test data in this repository will be used to test the pipeline from end-to-end. The associated parameters and settings to run the default tests for the pipeline can be found in [`test.config`](https://github.com/nf-core/methylong/blob/master/conf/test.config).

### Files

* `test_samplesheet.csv` - Sample information sheet required for the pipeline.
* `test_data/test_modbam/` - ONT and PacBio HiFi unaligned modified basecalled bam files (modbam) from a single individual Arabidopsis thaliana (Col-0) 2Mb region:  Chr1:5MB-7MB.

## Pod5 data

The pod5 data in this repository will be used to test the pipeline with the ONT basecalling step. The associated parameters and settings to run the default tests for the pipeline can be found in [`test_pod5.config`](https://github.com/nf-core/methylong/blob/master/conf/test_pod5.config).

### Files

* `test_samplesheet_pod5.csv` - Sample information sheet required for the pipeline.
* `v2.0.0/test_data/test_pod5/` - ONT pod5 files from HG002 5Mb region:  Chr1:5MB-7MB and 7MB-10MB.

## Unmodified BAM and m6A data

The unmodified BAM data in this repository will be used to test the pipeline with the Pacbio modcalling (include m6a-calling) step. The associated parameters and settings to run the default tests for the pipeline can be found in [`test_unmodified_bam.config`](https://github.com/nf-core/methylong/blob/master/conf/test_unmodified_bam.config).

### Files

* `test_samplesheet_unmodified_bam.csv` - Sample information sheet required for the pipeline
* `v2.0.0/test_data/test_bam/` - PacBio HiFi unaligned unmodified basecalled bam demo file from HG002.

## DMR population-scale data

The dmr data in this repository will be used to test the pipeline with the population-scale DMR step. The associated parameters and settings to run the pipeline can be found in [`test_dmr.config`](https://github.com/nf-core/methylong/blob/master/conf/test_dmr.config).

### Files

* `test_samplesheet_dmr.csv` - Sample information sheet required for the pipeline
* `v2.0.0/test_date/test_dmr/` - ONT pod5 files from HG002 and HG003 5Mb region:  Chr1:5MB-10MB.

## Full-sized data

The full-sized data in this repository will be used to test the pipeline from end-to-end. The associated parameters and settings to run the pipeline can be found in [`test_full.config`](https://github.com/nf-core/methylong/blob/master/conf/test_fullconfig).

### Files

* `full_test_samplesheet.csv` - Sample information sheet required for the pipeline

