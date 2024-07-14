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
for CHR in 21 22; do
    mkdir -p data/reference_genome/${CHR} data/panel/${CHR}
    # Download reference genome GRCh38
    wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr${CHR}.fa.gz -O data/reference_genome/${CHR}/hs38DH.chr${CHR}.fa.gz
    # Download phased panel
    wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz -O data/panel/${CHR}/1000GP.chr${CHR}.vcf.gz
    # Download phased panel index
    wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi -O data/panel/${CHR}/1000GP.chr${CHR}.vcf.gz.tbi
done

wget https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip -O data/reference_genome/GRCh38.map.zip
```

The affimetrix SNP array is also to be downloaded with

```bash
wget -c https://api.gdc.cancer.gov/data/9bd7cbce-80f9-449e-8007-ddc9b1e89dfb -O data/affi/snp6.txt.gz
gunzip data/affi/snp6.txt.gz
```

##### Merging

Datas can be separated by chromosome or not.
Hence we will aggregated both chr21 and chr22 into one.

```bash
# Merge the panel file
PANEL_21=./data/panel/21/1000GP.chr21.vcf.gz
PANEL_22=./data/panel/21_22/1000GP.chr21.vcf.gz
PANEL_21_22=./data/panel/21_22/1000GP.chr21_22.vcf.gz

bcftools concat -Oz -o ${PANEL_21_22} ${PANEL_21} ${PANEL_22}
bcftools index -f ${PANEL_21_22} --threads 4

# Merge the fasta file
FASTA_21=./data/reference_genome/21/hs38DH.chr21.fa
FASTA_22=./data/reference_genome/22/hs38DH.chr22.fa
FASTA_21_22=./data/reference_genome/21_22/hs38DH.chr21_22.fa
gunzip $FASTA_21.gz
gunzip $FASTA_22.gz
mkdir -p data/reference_genome/21_22
cat ${FASTA_21} ${FASTA_22} > ${FASTA_21_22}
```

##### Relationship analysis

To perform a representative imputation with the 1000 GP data we need to take out the individuals we want to impute from the panel as well as their relatives.
To do so we will perform with plink a relationship analysis with the `--genome` option. (This might take some minute to compute).

```bash
mkdir -p analysis
# Compute the relationship analysis with both chr21 and chr22
plink --vcf ${PANEL_21_22} \
    --mind 0.20 --geno 0.2 --maf 0.01 \
    --genome --out analysis/1000GP_chr21_22 

# Select the individuals to remove from the reference panel
Rscript --vanilla ./analysis/relationship_analysis.R \
    --input  ./analysis/1000GP_chr21_22.genome \
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
    analysis/selected_individuals.txt \
    analysis/ERR323_listing.txt
```

#### Preparation of the different reference files

1) Filter the region of interest from the fasta and the panel files
2) Filter the region of interest of the validation file gnomAD
3) Normalise the panel and filter out related individual to selected individuals
4) Select only the SNPS
5) Convert to TSV

```bash
for chr in 21 22 21_22; do
    echo $chr
    . get_panel_s.sh \
        data/panel/$chr/1000GP.chr$chr \
        data/reference_genome/$chr/hs38DH.chr$chr \
        region.lst
done
```

6) Extract the SNP position present in the SNP chip array
7) Get the map file for the reference genome

```bash
. get_map_snp.sh \
    data/reference_genome/ \
    GRCh38 \
    data/affi/snp6 \
    region.lst
```

#### Preparation and downsampling of the individual file validation and test file

1) Filter out the region of interest and format to BAM
2) Get the genotype likelihood based on the panel for the validation file and simulated file
3) Downsampling the individual data to 1X
4) Extract from the validation file the SNP position present in the SNP chip array

```bash
. get_ind_1x.sh \
    data/panel/21_22/1000GP.chr21_22.s.norel \
    data/reference_genome/21_22/hs38DH.chr21_22.fa \
    data/affi/snp6.s.map \
    region.lst
```

#### Impute with glimpse

```bash
. get_ind_imputed.sh \
    data/panel/21_22/1000GP.chr21_22.s.norel \
    data/reference_genome/21_22/hs38DH.chr21_22.fa \
    region.lst
```

### For Beagle

```bash
wget http://faculty.washington.edu/browning/beagle/test.22Jul22.46e.vcf.gz -O data/beagle/test.22Jul22.46e.vcf.gz

echo "*** Creating test files: ref.22Jul22.46e.vcf.gz target.22Jul22.46e.vcf.gz ***"
zcat data/beagle/test.22Jul22.46e.vcf.gz | cut -f1-190 | tr '/' '|' | gzip > data/beagle/ref.22Jul22.46e.vcf.gz
zcat data/beagle/test.22Jul22.46e.vcf.gz | cut -f1-9,191-200 | gzip > data/beagle/target.22Jul22.46e.vcf.gz
```

### For Stitch

The only additional file generated for this tool is the posfile, generated using bcftools from files already present in the modules branch. The resulting test file has been also put in the modules branch.

```bash
wget https://github.com/nf-core/test-datasets/raw/modules/data/genomics/homo_sapiens/genome/vcf/dbsnp_146.hg38.vcf.gz

bcftools query -i 'TYPE=="SNP" & N_ALT==1' -f '%CHROM\t%POS\t%REF\t%ALT' > dbsnp_146.hg38.biallelic_snps.tsv
```

## Files size

To get the size of all the files inside the git repository use the following command

```bash
git ls-files | xargs -r du -h | sort -h
```
