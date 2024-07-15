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

The aim of this dataset is to provided a minimal set allowing to impute multiple individuals across multiple chromosome and test the reliability of this imputation.
To do so we will use the data available from the [**1000 genome project**](http://ftp.1000genomes.ebi.ac.uk).
To reduce the size of the files only the region 16600000-16800000 of the chromosomes 21 and 22 of the human genome.

### Environment

To use the different script below you need bcftools, samtools, plink1.9 and tabix.
You can install everything with conda by using the following commands:

```bash
conda env create --name env_tools --file environment.yml
conda activate env_tools
```

### Initial data

#### Reference genome and Panel

##### Downloading

You first need to download the reference genome and phase panel for each chromosome.

```bash
for CHR in chr21 chr22; do
    mkdir -p hum_data/panel/${CHR}
    # Download phased panel
    wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz -O hum_data/panel/${CHR}/1000GP.${CHR}.vcf.gz
    # Download phased panel index
    wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi -O hum_data/panel/${CHR}/1000GP.${CHR}.vcf.gz.tbi
done

# Download reference genome GRCh38
mkdir -p hum_data/reference_genome/
wget -c -O- https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | gunzip | bgzip  > hum_data/reference_genome/GRCh38.fa.bgz

# Download the reference genome map
wget https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip -O hum_data/reference_genome/GRCh38.map.zip
```

The affimetrix SNP array is also to be downloaded with

```bash
wget -c https://api.gdc.cancer.gov/data/9bd7cbce-80f9-449e-8007-ddc9b1e89dfb -O hum_data/affi/snp6.txt.gz
```

##### Relationship analysis

To perform a representative imputation with the 1000 GP data we need to take out the individuals we want to impute from the panel as well as their relatives.
To do so we will perform with plink a relationship analysis with the `--genome` option. (This might take some minute to compute).

```bash
mkdir -p analysis
# Compute the relationship analysis on chr21 (takes a few minutes)
plink --vcf hum_data/panel/chr21/1000GP.chr21.vcf.gz \
    --mind 0.20 --geno 0.2 --maf 0.01 \
    --genome --out analysis/1000GP_chr21 

# Select the individuals to remove from the reference panel
Rscript --vanilla ./analysis/relationship_analysis.R \
    --input  ./analysis/1000GP_chr21.genome \
    --ind_sel ./ind_sel.txt \
    --output ./analysis
```

#### Individuals data

Individuals data, beware can be long to download
Individual NA12878 is from the Glimpse tutorial, NA19401 is one with only one relationship with a PI_HAT > 0.2
To find the corresponding folder use :

```bash
lftp -e "find;quit" ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323 > analysis/ERR323_listing.txt
```

We can now download the corresponding folder:

> [!WARNING]
> Each individuals `.cram` files are around 15Gb.

```bash
. download_ind.sh \
    ind_sel.lst \
    hum_data/individuals \
    analysis/ERR323_listing.txt \
    ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323
```

#### Preparation of the different reference files

1) Filter the region of interest from the fasta and the panel files
2) Normalise the panel and filter out related individual to selected individuals
3) Chunks the chromosomes
4) Create the sites file
5) Convert to haplegend format

```bash
. get_panel_s.sh \
    hum_data/panel \
    1000GP \
    hum_data/reference_genome/GRCh38 \
    region.lst \
    chr
```

6) Extract the SNP position present in the SNP chip array
7) Get the map file for the reference genome

```bash
. get_map_snp.sh \
    hum_data/reference_genome/ \
    GRCh38 \
    hum_data/affi/snp6 \
    region.lst
```

#### Preparation and downsampling of the individual file validation and test file

1) Filter out the region of interest and format to BAM
2) Downsampling the individual data to 1X

```bash
. get_ind_1x.sh \
    hum_data/individuals \
    hum_data/reference_genome/GRCh38.s.fa \
    ind_sel.lst \
    region.lst
```

3) Get the genotype likelihood based on the panel for the validation file and simulated file
4) Extract from the validation file the SNP position present in the SNP chip array

```bash
. get_ind_snp.sh \
    hum_data/individuals \
    hum_data/panel \
    1000GP \
    hum_data/reference_genome/GRCh38.s.fa.gz \
    hum_data/affi/snp6.s.map \
    ind_sel.lst \
    region.lst
```

#### Impute with Glimpse2

```bash
. get_ind_imputed.sh \
    hum_data/panel \
    1000GP \
    hum_data/individuals \
    hum_data/reference_genome/GRCh38.s.fa.gz \
    region.lst \
    ind_sel_1.lst
```

### For Beagle

```bash
wget http://faculty.washington.edu/browning/beagle/test.22Jul22.46e.vcf.gz -O data/beagle/test.22Jul22.46e.vcf.gz

echo "*** Creating test files: ref.22Jul22.46e.vcf.gz target.22Jul22.46e.vcf.gz ***"
zcat data/beagle/test.22Jul22.46e.vcf.gz | cut -f1-190 | tr '/' '|' | gzip > data/beagle/ref.22Jul22.46e.vcf.gz
zcat data/beagle/test.22Jul22.46e.vcf.gz | cut -f1-9,191-200 | gzip > data/beagle/target.22Jul22.46e.vcf.gz
```

## Files size

To get the size of all the files inside the git repository use the following command

```bash
git ls-files | xargs -r du -h | sort -h
```
