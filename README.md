# test-datasets: `raredisease`

This branch contains test data to be used for automated testing with the [nf-core/raredisease](https://github.com/nf-core/raredisease) pipeline.

## Content of this repository

`reference/`: background resources needed by tools of raredisease pipeline

`testdata/`: chr20 test resources

`testdata/grch38_gnomad_reformated_-r3.1.1-.vcf.gz`: Gnomad vcf file containing entries for the region chr20:90000-92000

`testdata/grch38_vcfanno_config_-v0.2-_chr20.toml`: TOML file for small test

`testdata/vcfanno_grch38_small_test.tar.gz`: the archived files of grch38_*.{vcf.gz, vcf.gz.tbi} for small test
