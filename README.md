# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

Test data to be used for automated testing with the nf-core pipelines

## Test data for nf-core/pacsomatic

This branch contains test data for the [nf-core/pacsomatic](https://github.com/nf-core/pacsomatic) pipeline.

## Introduction

The nf-core/pacsomatic pipeline is designed for somatic variant calling using PacBio long-read sequencing data. This test dataset provides a minimal set of files to validate the pipeline's functionality during automated testing and development.

## Test Dataset Composition

### Source Data

Original data: [GIAB HG008 Cancer Genome in a Bottle](https://www.nist.gov/programs-projects/cancer-genome-bottle)  
Download location: [GIAB FTP - PacBio Revio 20240125](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/PacBio_Revio_20240125/)

### Sample Information

| Patient       | Sample ID | Status | Data Type                          |
| ------------- | --------- | ------ | ---------------------------------- |
| Patient_HG008 | DS_MT_T   | Tumor  | Downsampled mitochondrial BAM      |
| Patient_HG008 | DS_MT_N   | Normal | Downsampled mitochondrial BAM      |

### Rationale

The mitochondrial genome (16.5 kb) was selected for fast CI/CD testing while maintaining sufficient complexity for somatic variant calling validation. Mitochondrial heteroplasmy provides somatic-like variant patterns suitable for pipeline testing.

### Files Included

#### Sequencing Data (`testdata/`)

Downsampled BAM files containing mitochondrial reads generated from GIAB HG008 samples. 

**Generation procedure (example for normal sample):**

1. Downsample to ~10K unaligned reads:
   ```bash
   samtools view -c -f 4 ${HG008_Normal_hifibam}  # Get total reads
   samtools view -b -s 0.00015 ${HG008_Normal_hifibam} > HG008-N_Downsample_10K.bam
   ```

2. Align to hg38 using pbmm2:
   ```bash
   pbmm2 align hg38.fa HG008-N_Downsample_10K.bam Patient_HG008_Downsample_10K_N.align.bam
   ```

3. Extract mitochondrial reads:
   ```bash
   samtools view -b Patient_HG008_Downsample_10K_N.align.bam MT | \
     samtools sort -@ 16 -o HG008_Downsample_MT_normal.align.bam
   ```

4. Revert to unaligned BAM:
   ```bash
   java -jar picard.jar RevertSam \
     INPUT=HG008_Downsample_MT_normal.align.bam \
     OUTPUT=HG008_Downsample_MT_normal.bam
   ```

**Tools used:** samtools/1.17, picard/2.9.4, pbmm2 (smrttools/25.1)

The tumor sample was generated using the same procedure.

#### Reference Files (`reference/`)

- `MT.fa`: Human mitochondrial genome reference sequence ([NCBI J01415.2](https://www.ncbi.nlm.nih.gov/nuccore/J01415.2), 16,569 bp)
- `MT_target.bed`: Manually created BED file defining target regions covering the mitochondrial genome (MT:0-3106, MT:3107-16569)

### Samplesheet Format

The `samplesheet.csv` follows the nf-core/pacsomatic input schema:

```csv
patient,sample,status,bam,pbi
Patient_HG008,DS_MT_T,1,https://raw.githubusercontent.com/nf-core/test-datasets/pacsomatic/testdata/HG008_Downsample_MT_tumor.bam,
Patient_HG008,DS_MT_N,0,https://raw.githubusercontent.com/nf-core/test-datasets/pacsomatic/testdata/HG008_Downsample_MT_normal.bam,
```

**Column descriptions:**
- `patient`: Patient identifier (must be unique per patient)
- `sample`: Sample identifier (must be unique)
- `status`: Sample status (0 = normal, 1 = tumor)
- `bam`: Path or URL to the BAM file
- `pbi`: Optional path to PacBio index file (.pbi)

## Documentation

nf-core/test-datasets comes with documentation in the `docs/` directory:

1. [Add a new test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/ADD_NEW_DATA.md)
2. [Use an existing test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/USE_EXISTING_DATA.md)

## Downloading test data

Due the large number of large files in this repository for each pipeline, we highly recommend cloning only the branches you would use.

```bash
git clone <url> --single-branch --branch pacsomatic
```

To subsequently clone other branches[^1]

```bash
git remote set-branches --add origin [remote-branch]
git fetch
```

## Usage with nf-core/pacsomatic

To run the pipeline with this test dataset:

```bash
nextflow run nf-core/pacsomatic \
  -profile test,docker \
  --outdir results
```

Or with a local clone:

```bash
nextflow run nf-core/pacsomatic \
  -profile docker \
  --input samplesheet.csv \
  --fasta reference/MT.fa \
  --target_bed reference/MT_target.bed \
  --outdir results
```

## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).

[^1]: From [stackoverflow](https://stackoverflow.com/a/60846265/11502856)
