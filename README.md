# test-datasets `eager`
Test data to be used for automated testing with the nf-core pipelines

## Content of this repository

For both `reference` and `testdata`, there is one directory for each organism type. 

In `references` this will contain fasta, possible indices and auxiliary files (such as annotation files). 

In `testdata`, there will be one directory for file type, e.g. `fastq`, `bam`, `vcf`. 

### Reference Genomes (`reference/`)

#### Mammoth 

`Mammoth_MT_Krause.fasta`: Reference genome for mammoth data, plus bwa/samtools/picard indices of this file.
`Mammoth_MT_Krause.gff3`: GFF file with feature annotations for the above FASTA file.

#### Human

`1240K.pos.list_hs37d5.0based.bed.gz`: a bed file containing positions for the '1240k' SNP capture array (file by [Alex Peltzer](https://github.com/apeltzer) and [Stephen Clayton](https://github.com/sc13-bioinf), originally defined in [Mathieson et al. 2015 _Nature_](https://www.nature.com/articles/nature16152))

`1240K.pos.list_GH19.0based.bed.gz`: a bed file containing positions for the '1240k' SNP capture array for the HG19 reference genome ported from the above by [Thiseas C. Lamnidis](https://github.com/TCLamnidis), originally defined in [Mathieson et al. 2015 _Nature_](https://www.nature.com/articles/nature16152))

### Sequencing Data (`testdata/`)

#### Mammoth

TSV input versions of FASTQ and BAM files exists for all files below.

**fastq**

This contains two paired end samples which are the default test samples.

`JK2782_*`: HiSeq MT captured library with no UDG treatment (~10K reads after merging, [Fellows Yates et al. 2017 _Sci. Rep._](https://doi.org/10.1038/s41598-017-17723-1))
`JK2802_*`: HiSeq MT captured library with no UDG treatment (~10K reads after merging, [Fellows Yates et al. 2017 _Sci. Rep._](https://doi.org/10.1038/s41598-017-17723-1))

**bam**

Already clipped, merged and mapped BAM files of the two paired-end FASTQ files as described above in `fastq`.

**vcf**

`JK2772_*`: An additional VCF for another Mammoth MT capture library that was previously processed in nf-core/eager and genotyped via UnifiedGenotyper (~10K reads after merging, [Fellows Yates 2017 et al._Sci. Rep._](https://doi.org/10.1038/s41598-017-17723-1))

**maltextract**

`MaltExtract_list.txt`: A list of taxa in the MALT database (see [below](#malt)), for running MaltExtract.

#### Human

TSV input versions of BAM files exists for all files below.

**bam**

`JK2067_*`: HiSeq 1240k captured UDG-half single-end (~10K reads after clipping, [Lamnidis et al. 2018 _Nat. Comms._](https://doi.org/10.1038/s41467-018-07483-5))

#### Benchmarking

There are three test TSV input files that can be used for larger, more 'realistic' shotgun testing with full sized data. This currently covers three main contexts.

`benchmarking_human.tsv`: Ancient Fish from [Star et al. 2017](https://doi.org/10.1073/pnas.1710186114)
`benchmarking_vikingfish.tsv`: Ancient Humans from [Gamba et al. 2014](https://doi.org/10.1073/10.1038/ncomms6257)
`benchmarking_pathogenscreening.tsv`: Ancient Pathogen from [Andrades Valtueña et al. 2017](https://doi.org/10.1016/j.cub.2017.10.025)

A further one for internal large-scale AWS testing is:

`human_stresstest.tsv`: 137 ancient human genomes from [de Barros Damgaard et al. 2018](https://doi.org/10.1038/s41586-018-0094-2)

### Databases

#### MALT

Contains database for MALT, built with MALT v041, containing
  - A mammoth mitochondrial genome (Mammoth_MT_Krause above) 
  - A human mitochondrial genome (NCBI Accession: NC_012920.1)
  - A boa constrictor mitochondrial genome (NCBI Accession: AB177354.1)


#### Kraken

Contains database for Kraken2, built with Kraken 2.0.8-beta, containing
  - A mammoth mitochondrial genome (NCBI Accession: NC_007596.2) 
  - A human mitochondrial genome (NCBI Accession: NC_012920.1)
  - A boa constrictor mitochondrial genome (NCBI Accession: AB177354.1)
