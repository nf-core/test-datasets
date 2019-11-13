# test-datasets `eager`
Test data to be used for automated testing with the nf-core pipelines

## Content of this repository

For both `reference` and `testdata`, there is one directory for each organism type. 

In `references` this will contain fasta, possible indices and auxiliary files (such as annotation files). 

In `testdata`, there will be one directory for file type, e.g. `fastq`, `bam`, `vcf`. 

### Reference genome(s)

#### Mammoth 

`Mammoth_MT_Krause.fasta`: Reference genome for mammoth data, plus bwa/samtools/picard indices of this file.
`Mammoth_MT_Krause.gff3`: GFF file with feature annotations for the above FASTA file.

#### Human

`1240K.pos.list_hs37d5.0based.bed.gz`: a bed file containing positions for the '1240k' SNP capture array (file by [Alex Peltzer](https://github.com/apeltzer) and [Stephen Clayton](https://github.com/sc13-bioinf), originally defined in [Mathieson et al 2015 _Nature_](https://www.nature.com/articles/nature16152))

### testdata

#### Mammoth

**fastq**

This contains two paired end samples which are the default test samples.

`JK2782_*`: HiSeq MT captured library with no UDG treatment (~10K reads after merging, [Fellows Yates et al. 2017 _Sci. Rep._](https://doi.org/10.1038/s41598-017-17723-1))
`JK2802_*`: HiSeq MT captured library with no UDG treatment (~10K reads after merging, [Fellows Yates et al. 2017 _Sci. Rep._](https://doi.org/10.1038/s41598-017-17723-1))

**bam**

Already clipped, merged and mapped BAM files of the two paired-end FASTQ files as described above in `fastq`.

**vcf**

`JK2772_*`: An additional VCF for another Mammoth MT capture library that was previously processed in nf-core/eager and genotyped via UnifiedGenotyper (~10K reads after merging, [Fellows Yates (2017) Sci. Rep](https://doi.org/10.1038/s41598-017-17723-1))

#### Human

**bam**

`JK2067_*`: HiSeq 1240k captured UDG-half single-end (~10K reads after clipping, [Lamnidis et al. 2018 _Nat. Comms._](https://doi.org/10.1038/s41467-018-07483-5))
