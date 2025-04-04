# Spaceranger test data

Based on example data from 10x Genomics website

## Human ovarian cancer 1, FFPE direct placement, v1 chemistry

Source: https://www.10xgenomics.com/resources/datasets/human-ovarian-cancer-1-standard

The example data has been downloaded and subsampled using the following bash script

```bash
#!/bin/bash
# from https://www.10xgenomics.com/resources/datasets/human-ovarian-cancer-1-standard
DIR=human-ovarian-cancer-1-standard_v1_ffpe
mkdir -p $DIR && cd $DIR

# Input Files
curl -O https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Ovarian_Cancer/Visium_FFPE_Human_Ovarian_Cancer_fastqs.tar
curl -O https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Ovarian_Cancer/Visium_FFPE_Human_Ovarian_Cancer_image.jpg
curl -O https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Ovarian_Cancer/Visium_FFPE_Human_Ovarian_Cancer_probe_set.csv

# Extract
tar xvf Visium_FFPE_Human_Ovarian_Cancer_fastqs.tar

# Create subsampled dataset with ImageMagick
# https://imagemagick.org/index.php
mkdir subsampled
convert Visium_FFPE_Human_Ovarian_Cancer_image.jpg -resize 1500x1500 subsampled/Visium_FFPE_Human_Ovarian_Cancer_image.jpg
for f in Visium_FFPE_Human_Ovarian_Cancer_fastqs/*L001*R*; do; gzip -cdf $f | head -n 40000 | gzip -c > subsampled/$(basename $f); done
```

## Human brain cancer, FFPE cytassist, v2 chemistry

Source: https://www.10xgenomics.com/resources/datasets/human-brain-cancer-11-mm-capture-area-ffpe-2-standard

The example data has been downloaded and subsampled using the following bash script

```bash
DIR=human-brain-cancer-11-mm-capture-area-ffpe-2-standard_v2_ffpe_cytassist
mkdir -p $DIR && cd $DIR

# Input Files
curl -O https://s3-us-west-2.amazonaws.com/10x.files/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Glioblastoma/CytAssist_11mm_FFPE_Human_Glioblastoma_fastqs.tar
curl -O https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Glioblastoma/CytAssist_11mm_FFPE_Human_Glioblastoma_image.tif
curl -O https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Glioblastoma/CytAssist_11mm_FFPE_Human_Glioblastoma_probe_set.csv
curl -O https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Glioblastoma/CytAssist_11mm_FFPE_Human_Glioblastoma_tissue_image.tif

tar xvf CytAssist_11mm_FFPE_Human_Glioblastoma_fastqs

mkdir subsampled
# cytassist only takes original image... 24M is still manageable, even for a test dataset
cp CytAssist_11mm_FFPE_Human_Glioblastoma_image.tif subsampled/
for f in CytAssist_11mm_FFPE_Human_Glioblastoma_fastqs/*S1*L00{1,2}*R*; do; gzip -cdf $f | head -n 40000 | gzip -c > subsampled/$(basename $f); done
```

## Human lung cancer, FFPE cytassist, HD only

Source: https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-human-lung-cancer-post-xenium-expt

The example data has been downloaded and subsampled using the following bash script:

```bash
DIR=visium-hd-human-lung-cancer-hd-ffpe-2-standard_v2_ffpe_cytassist
mkdir -p $DIR && cd $DIR

# Input Files
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer_HD_Only_Experiment2/Visium_HD_Human_Lung_Cancer_HD_Only_Experiment2_image.tif
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer_HD_Only_Experiment2/Visium_HD_Human_Lung_Cancer_HD_Only_Experiment2_tissue_image.btf
curl -O https://s3-us-west-2.amazonaws.com/10x.files/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer_HD_Only_Experiment2/Visium_HD_Human_Lung_Cancer_HD_Only_Experiment2_fastqs.tar
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Lung_Cancer_HD_Only_Experiment2/Visium_HD_Human_Lung_Cancer_HD_Only_Experiment2_probe_set.csv
url -O https://cf.10xgenomics.com/supp/spatial-exp/refdata-gex-GRCh38-2020-A.tar.gz

# Extract
tar xvf Visium_HD_Human_Lung_Cancer_HD_Only_Experiment2_fastqs.tar

# Create subsampled dataset with ImageMagick
# https://imagemagick.org/index.php
mkdir subsampled
convert Visium_HD_Human_Lung_Cancer_HD_Only_Experiment2_tissue_image.btf -compress JPEG -quality 1 subsampled/Visium_HD_Human_Lung_Cancer_HD_Only_Experiment2_tissue_image.btf
cp Visium_HD_Human_Lung_Cancer_HD_Only_Experiment2_image.tif subsampled/

# Untar the reference genome
tar xvfz refdata-gex-GRCh38-2020-A.tar.gz

# Run spaceranger count to get bam file
spaceranger count \
          --id "Visium_HD_Human_Lung_Cancer_HD_Only_D" \
          --sample "Visium_HD_Human_Lung_Cancer_HD_Only_D" \
          --fastqs="./Visium_HD_Human_Lung_Cancer_HD_Only_Experiment2_fastqs" \
          --slide "H1-84QJZFR" \
          --area "D1" \
          --transcriptome="./refdata-gex-GRCh38-2020-A" \
          --probe-set="Visium_HD_Human_Lung_Cancer_HD_Only_Experiment2_probe_set.csv" \
          --image="Visium_HD_Human_Lung_Cancer_HD_Only_Experiment2_tissue_image.btf" \
          --cytaimage="Visium_HD_Human_Lung_Cancer_HD_Only_Experiment2_image.tif" \
          --create-bam true \
          --output-dir "./outs"

# Create filtered_barcodes.csv from loupe output for cropping subregion.

# Subsample the bam file to get a smaller test dataset, based on a set of barcodes
subset-bam --bam ./outs/possorted_genome_bam.bam --cell-barcodes ./filtered_barcodes.csv --out-bam ./outs/cropped_bam.bam
bamtofastq ./outs/cropped_bam.bam ./subsampled
```