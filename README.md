# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the nf-core pipelines

> ⚠️ **Do not merge your test data to `master`! Each pipeline has a dedicated branch (and a special one for modules)**

## Introduction

nf-core is a collection of high quality Nextflow pipelines. This repository contains various files for CI and unit testing of nf-core pipelines and infrastructure.

The principle for nf-core test data is as small as possible, as large as necessary. Please see the [guidelines](https://nf-co.re/docs/contributing/test_data_guidelines) for more detailed information. Always ask for guidance on the [nf-core slack](https://nf-co.re/join) before adding new test data.

## Documentation

nf-core/test-datasets comes with documentation in the `docs/` directory:

01. [Add a new  test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/ADD_NEW_DATA.md)
02. [Use an existing test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/USE_EXISTING_DATA.md)

## Downloading test data

Due the large number of large files in this repository for each pipeline, we highly recommend cloning only the branches you would use.

```bash
git clone <url> --single-branch --branch <pipeline/modules/branch_name>
```

To subsequently clone other branches[^1]

```bash
git remote set-branches --add origin [remote-branch]
git fetch
```

## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).

[^1]: From [stackoverflow](https://stackoverflow.com/a/60846265/11502856)

## How data has been acquired

### For GLIMPSE

#### Initial data

For CHR 21

```
mkdir -p data/reference_genome
wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.filtered.shapeit2-duohmm-phased.vcf.gz -O data/panel/panel_2020-08-05_chr21.phased.vcf.gz
wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.filtered.shapeit2-duohmm-phased.vcf.gz.tbi -O data/panel/panel_2020-08-05_chr21.phased.vcf.gz.tbi
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz -O data/reference_genome/hs38DH.chr21.fa.gz
```

For CHR 22

```
wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz -O data/panel/panel_2020-08-05_chr22.phased.vcf.gz
wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz.tbi -O data/panel/panel_2020-08-05_chr22.phased.vcf.gz.tbi
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz -O data/reference_genome/hs38DH.chr22.fa.gz
```

Individuals data, beware can be long to download
Individual NA12878 is from the Glimpse tutorial, NA19401 is one with only one relationship with a PI_HAT > 0.2
To find the corresponding folder use :

```
lftp -e "find;quit" ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323 > listing.txt
```

We can now download the corresponding folder:

```
mkdir -p data/NA12878 data/NA19401
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram -O data/NA12878/NA12878.final.cram
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram.crai -O data/NA12878/NA12878.final.cram.crai
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239749/NA19401.final.cram -O data/NA19401/NA19401.final.cram
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239749/NA19401.final.cram.crai -O data/NA19401/NA19401.final.cram.crai
```

SNP array value

```
wget -c https://api.gdc.cancer.gov/data/9bd7cbce-80f9-449e-8007-ddc9b1e89dfb -O data/affi/snp6.txt.gz
gunzip data/affi/snp6.txt.gz
```

#### Environment

To use the different script below you need bcftools, samtools and tabix.
You can install everything with conda by using the following commands:

```
conda env create --name env_tools --file environment.yml
conda activate env_tools
```

#### Preparation of the different panel files

1) Filter the region of interest of the panel file
2) Filter the region of interest of the validation file gnomAD
3) Normalise the panel and filter out related individual to NA12878
4) Select only the SNPS
5) Convert to TSV

```
. get_panel_s.sh
```

#### Preparation and downsampling of the individual file validation and test file

1) Filter out the region of interest and format to BAM
2) Get the genotype likelihood based on the panel for the validation file
3) Downsampling the individual data to 1X

```
. get_ind_1x
```

#### Compute the genotype likelihood for the individual data 

1) Compute genotype likelihood based on the panel

```
. get_ind_gl.sh
```

### For Beagle

```
wget http://faculty.washington.edu/browning/beagle/test.22Jul22.46e.vcf.gz

echo "*** Creating test files: ref.22Jul22.46e.vcf.gz target.22Jul22.46e.vcf.gz ***"
zcat test.22Jul22.46e.vcf.gz | cut -f1-190 | tr '/' '|' | gzip > ref.22Jul22.46e.vcf.gz
zcat test.22Jul22.46e.vcf.gz | cut -f1-9,191-200 | gzip > target.22Jul22.46e.vcf.gz

wget https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
```

### For Stitch

The only additional file generated for this tool is the posfile, generated using bcftools from files already present in the modules branch. The resulting test file has been also put in the modules branch.
```
wget https://github.com/nf-core/test-datasets/raw/modules/data/genomics/homo_sapiens/genome/vcf/dbsnp_146.hg38.vcf.gz

bcftools query -i 'TYPE=="SNP" & N_ALT==1' -f '%CHROM\t%POS\t%REF\t%ALT' > dbsnp_146.hg38.biallelic_snps.tsv
```

## Files size

To get the size of all the files inside the git repository use the following command
```
git ls-files | xargs -r du -h | sort -h
```
