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

NA11995 HapMap CEPH HLA-fished reads (SRA: SRR766010 — sourced from this branch's
`fastq/NA11995_SRR766010_*_fished.fastq.gz`) re-paired with seqkit and aligned to
GRCh38 with `bwa mem`. Used by hlatyping's BAM-input test profiles (OptiType,
HLA*LA, SpecHLA).

Reference: 1000 Genomes `GRCh38_full_analysis_set_plus_decoy_hla.fa` ([FTP](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa)).

```bash
# Re-pair the unsynchronised fished FASTQs (R1 and R2 have different read sets)
seqkit pair -1 NA11995_SRR766010_1_fished.fastq.gz \
            -2 NA11995_SRR766010_2_fished.fastq.gz \
            -O paired/

# Index reference (one-time, ~75 min)
bwa index GRCh38_full_analysis_set_plus_decoy_hla.fa

# Align, sort, index
bwa mem -R '@RG\tID:NA11995\tSM:NA11995' \
        GRCh38_full_analysis_set_plus_decoy_hla.fa \
        paired/NA11995_R1.fastq.gz paired/NA11995_R2.fastq.gz \
  | samtools sort -o NA11995_GRCh38.bam -
samtools index NA11995_GRCh38.bam
```

## Support

For further information or help, don't hesitate to get in touch on our [Gitter channel](https://gitter.im/nf-core/Lobby)
