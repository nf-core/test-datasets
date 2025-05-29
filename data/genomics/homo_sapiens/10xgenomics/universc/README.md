# [universc](https://github.com/minoda-lab/universc)

This tool used an old version of cellranger for its reference.
Therefore you need to use the `nf-core/universc:1.2.5.1` container with it's
embedded cellranger to generate the correct reference layout to avoid errors.

The following script has been used

```bash
cellranger \
    mkgtf \
    data/genomics/homo_sapiens/genome/genome.gtf \
    data/genomics/homo_sapiens/10xgenomics/universc/genome.filtered.gtf

cellranger \
    mkref \
    --genome=homo_sapiens_chr22_reference \
    --fasta=data/genomics/homo_sapiens/genome/genome.fasta \
    --genes=data/genomics/homo_sapiens/10xgenomics/universc/genome.filtered.gtf

tar -czvf homo_sapiens_chr22_reference.tar.gz homo_sapiens_chr22_reference

rm data/genomics/homo_sapiens/10xgenomics/universc/genome.filtered.gtf
```
