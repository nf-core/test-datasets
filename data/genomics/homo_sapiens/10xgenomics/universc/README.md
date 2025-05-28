# [universc](https://github.com/minoda-lab/universc)

This tool used an old version of cellranger for its reference.
Therefore you need to use the `nf-core/universc:1.2.5.1` container with it's
embedded cellranger to generate the correct reference layout to avoid errors.

The following script has been used

```
cellranger \
    mkgtf \
    genome.gtf \
    genome.filtered.gtf

cellranger \
    mkref \
    --genome=homo_sapiens_chr22_reference \
    --fasta=genome.fasta \
    --genes=genome.filtered.gtf

tar -czvf homo_sapiens_chr22_reference.tar.gz homo_sapiens_chr22_reference
```
