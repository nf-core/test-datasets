# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

Test data to be used for automated testing with the nf-core pipelines

## Introduction

nf-core is a collection of high quality Nextflow pipelines.

## Documentation

nf-core/test-datasets comes with documentation in the `docs/` directory:

01. [Add a new  test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/ADD_NEW_DATA.md)
02. [Use an existing test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/USE_EXISTING_DATA.md)

## Data generation

### STARsolo / AlevinQC Testdata

Please ask [Olga Botvinnik](https://github.com/olgabot) for details on how this data was generated and subsetted.

### Kallisto/Bustools Testdata

The [reference/kallisto] and [testdata/kallisto] folders hold testing data that was subsetted to be able to utilize the data on automated continous integration services due to memory and time restrictions on these services. The data used here refers to [this howto article](https://www.kallistobus.tools/getting_started_explained.html). The files were subsetted utilizing these commands:

```bash
zcat SRR8599150_S1_L001_R1_001.fastq.gz |head -n 5000 > SRR8599150_S1_L001_R1_001.sub5000.fastq
zcat SRR8599150_S1_L001_R2_001.fastq.gz |head -n 5000 > SRR8599150_S1_L001_R2_001.sub5000.fastq
zcat Mus_musculus.GRCm38.cdna.all.fa.gz | sed -n '433032,517910 p'
zcat Mus_musculus.GRCm38.96.gtf.gz | grep -e '^#' -e '^19' > Mus_musculus.GRCm38.96.chr19.gtf
```

The GTF file contains annotation for more than just the chr19 data but has large portions of exons on chr19, so gives somewhat meaningful results. The cdna file was evaluated manually to determine an appropriate range with chr19 entries for testing.

## Support

For further information or help, don't hesitate to get in touch on our [Slack](https://nfcore.slack.com) or [Click here for an invite](https://nf-core-invite.herokuapp.com/)
