# test-datasets `eager2`
Test data to be used for automated testing with the nf-core pipelines

## Content of this repository

For both `reference` and `testdata`, there is one directory for each organism type. 

In `references` this will contain fasta, possible indices and auxiliary files (such as annotation files). 

In `testdata`, there will be one directory for file type, e.g. `fastq`, `bam`, `vcf`. 

### Reference genome(s)

#### Mammoth 

`Mammoth_MT_Krause.fasta`: Reference genome for mammoth data.

### testdata

#### Mammoth

**fastq**

This contains two paired end samples which are the default test samples.

`JK2782_*`: HiSeq MT captured library with no UDG treatment (~10K reads after merging, Fellows-Yates (2017) Sci. Rep)
`JK2802_*`: HiSeq MT captured library with no UDG treatment (~10K reads after merging, Fellows-Yates (2017) Sci. Rep)

**bam**

Already clipped, merged and mapped BAM files of the two paired-end FASTQ files as described above in [fastq](#fastq).

**vcf**

`JK2772_*`: An additional VCF for another Mammoth MT capture library that was previously processed in nf-core/eager and genotyped via UnifiedGenotyper (~10K reads after merging, Fellows-Yates (2017) Sci. Rep)

#### Human

**bam**

`JK2067_*`: HiSeq 1240k captured UDG-half single-end (~10K reads after clipping, Lamnidis (2018) Nat. Comms.)
