# test-datasets: `pixelator`

This branch contains test data to be used for automated testing with the [nf-core/pixelator](https://github.com/nf-core/pixelator) pipeline.

## Content of this repository

- `testdata/` : subsampled fastq files
- `panels/`: panel files to test passing custom panels

- `samplesheet/samplesheet.csv`: Experiment design file for minimal test dataset.
- `samplesheet/samplesheet_v2.csv`: Experiment design file for minimal test dataset using the v2 panel.
- `samplesheet/samplesheet_full.csv`: Experiment design file for full test dataset.

## Test datasets origin

Molecular pixelation data retrieved from public datasets provided by Pixelgen Technologies AB.

- https://software.pixelgen.com/datasets/1k-human-pbmcs-v1.0-immunology-I
- https://software.pixelgen.com/datasets/uropod-t-cells-v1.0-immunology-I

Immunology II datasets are taken from:

s3://pixelgen-technologies-datasets/mpx-datasets/scsp/2.0/5-donors-pbmcs-v2.0/


### Sampling procedure


Data was subsampled using `seqkit sample` to 200k or 300k paired-end reads.
The default random seed was used.

```bash
seqkit sample -2 -n 300000 Uropod_control_R1_001.fastq.gz -o uropod_300k_control_R1_001.fastq.gz
seqkit sample -2 -n 300000 Uropod_control_R2_001.fastq.gz -o uropod_300k_control_R2_001.fastq.gz
seqkit sample -2 -n 200000 Sample01_human_pbmcs_unstimulated_200k_R1_001.fastq.gz -o Sample01_human_pbmcs_unstimulated_200k_R1_001.fastq.gz
seqkit sample -2 -n 200000 Sample01_human_pbmcs_unstimulated_200k_R2_001.fastq.gz -o Sample01_human_pbmcs_unstimulated_200k_R2_001.fastq.gz
seqkit sample -2 -n 500000 Sample01_PBMC_1085_r1_R1_001.fastq.gz -o v2/Sample01_PBMC_1085_500k_r1_R1_001.fastq.gz
seqkit sample -2 -n 500000 Sample01_PBMC_1085_r1_R2_001.fastq.gz -o v2/Sample01_PBMC_1085_500k_r1_R2_001.fastq.gz
```

To test input sample concatenation the uropod_control sample was split using `seqkit split2`

```bash
seqkit split2 -p 2 --read1 uropod_control_300k_R1_001.fastq.gz --read2 uropod_control_300k_R2_001.fastq.gz
```

Seqkit is available on bioconda.