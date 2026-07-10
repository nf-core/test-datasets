# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the nf-core pipelines

## Introduction

nf-core is a collection of high quality Nextflow pipelines.

## Documentation
nf-core/test-datasets comes with documentation in the `docs/` directory:

01. [Add a new  test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/ADD_NEW_DATA.md)
02. [Use an existing test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/USE_EXISTING_DATA.md)

## hlatyping test data

### `bam/NA11995_GRCh38.bam`

NA11995 reads from `fastq/NA11995_SRR766010_*_fished.fastq.gz`, re-paired and aligned to GRCh38.

```bash
seqkit pair -1 NA11995_SRR766010_1_fished.fastq.gz \
            -2 NA11995_SRR766010_2_fished.fastq.gz -O paired/

# Reference: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
bwa index GRCh38_full_analysis_set_plus_decoy_hla.fa
bwa mem -R '@RG\tID:NA11995\tSM:NA11995' \
        GRCh38_full_analysis_set_plus_decoy_hla.fa \
        paired/NA11995_R1.fastq.gz paired/NA11995_R2.fastq.gz \
  | samtools sort -o NA11995_GRCh38.bam -
samtools index NA11995_GRCh38.bam
```

### `samplesheets/samplesheet_full.csv`

Full-size input for `-profile test_full`. NA12878 and GM12878 are the same individual, so DNA and RNA calls can be cross-checked.

| Sample | Data | Source |
| ------ | ---- | ------ |
| `NA12878_WES` | Agilent WES FASTQ | [Zenodo 6513789](https://zenodo.org/record/6513789), shared with nf-core/sarek's full test |
| `GM12878_RNA` | RNA-seq FASTQ (SRX1603629) | `s3://ngi-igenomes/test-data/rnaseq/`, shared with nf-core/rnaseq's full test |
| `HepG2_PEPTIDE` | `peptide/HepG2_A.tsv` | Peptide input for immunotype |

## Support

For further information or help, don't hesitate to get in touch on our [Gitter channel](https://gitter.im/nf-core/Lobby)
