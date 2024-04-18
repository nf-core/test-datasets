# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

This branch contains test datasets as well as various files for CI and unit
testing for the [`nf-core/spatialtranscriptomics`](https://github.com/nf-core/spatialtranscriptomics)
pipeline. The test datasets comes from two separate samples obtained from the
official 10x Genomics website, which have subsequently been sub-sampled for
testing purposes.

## Contents

The test data includes the following:

- Raw FASTQ files
- Image files
- Probe sets
- Samplesheets (both with and without Space Ranger pre-processing)
- Processed Space Ranger data for downstream testing
- Compressed archives for downloading during CI testing

## Origins

### Human ovarian cancer 1, FFPE direct placement, v1 chemistry

Source: https://www.10xgenomics.com/resources/datasets/human-ovarian-cancer-1-standard

```bash
#!/bin/bash
DIR=human-ovarian-cancer-1-standard_v1_ffpe
mkdir -p $DIR && cd $DIR

# Input Files
curl -O https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Ovarian_Cancer/Visium_FFPE_Human_Ovarian_Cancer_fastqs.tar
curl -O https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Ovarian_Cancer/Visium_FFPE_Human_Ovarian_Cancer_image.jpg
curl -O https://cf.10xgenomics.com/samples/spatial-exp/1.3.0/Visium_FFPE_Human_Ovarian_Cancer/Visium_FFPE_Human_Ovarian_Cancer_probe_set.csv

# Extract
tar xvf Visium_FFPE_Human_Ovarian_Cancer_fastqs.tar

# Create sub-sampled dataset with ImageMagick
# https://imagemagick.org/index.php
mkdir sub-sampled
convert Visium_FFPE_Human_Ovarian_Cancer_image.jpg -resize 1500x1500 sub-sampled/Visium_FFPE_Human_Ovarian_Cancer_image.jpg
for f in Visium_FFPE_Human_Ovarian_Cancer_fastqs/*L001*R*; do; gzip -cdf $f | head -n 40000 | gzip -c > sub-sampled/$(basename $f); done
```

### Human brain cancer, FFPE cytassist, v2 chemistry

Source: https://www.10xgenomics.com/resources/datasets/human-brain-cancer-11-mm-capture-area-ffpe-2-standard

```bash
DIR=human-brain-cancer-11-mm-capture-area-ffpe-2-standard_v2_ffpe_cytassist
mkdir -p $DIR && cd $DIR

# Input Files
curl -O https://s3-us-west-2.amazonaws.com/10x.files/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Glioblastoma/CytAssist_11mm_FFPE_Human_Glioblastoma_fastqs.tar
curl -O https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Glioblastoma/CytAssist_11mm_FFPE_Human_Glioblastoma_image.tif
curl -O https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Glioblastoma/CytAssist_11mm_FFPE_Human_Glioblastoma_probe_set.csv
curl -O https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Glioblastoma/CytAssist_11mm_FFPE_Human_Glioblastoma_tissue_image.tif

tar xvf CytAssist_11mm_FFPE_Human_Glioblastoma_fastqs

mkdir sub-sampled
# cytassist only takes original image... 24M is still manageable, even for a test dataset
cp CytAssist_11mm_FFPE_Human_Glioblastoma_image.tif sub-sampled/
for f in CytAssist_11mm_FFPE_Human_Glioblastoma_fastqs/*S1*L00{1,2}*R*; do; gzip -cdf $f | head -n 40000 | gzip -c > sub-sampled/$(basename $f); done
```

## Support

For further information or help, don't hesitate to get in touch on our
[Slack organisation](https://nf-co.re/join/slack) (a tool for instant
messaging).
