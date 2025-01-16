# Test data set provenance

This test data set was generated as test data for the NBISweden/Earth-Biogenome-Project-pilot 
pipeline, located at
https://github.com/NBISweden/Earth-Biogenome-Project-pilot/tree/main/tests/data/tiny

## Authors

- Guilherme Dias
- Martin Pippel

## Provenance

- species: _Drosophila melanogaster_
- PacBio HiFi data was downloaded from NBCI´s SRA [SRR10238607](https://www.ncbi.nlm.nih.gov/sra/SRR10238607)
- Illumina HiC data was downloaded from NBCI´s SRA [SRR10512944](https://www.ncbi.nlm.nih.gov/sra/SRR10512944)

- PacBio HiFi data were randomly subsampled to 15X read coverage and assembled with [hifiasm version 0.19.8](https://github.com/chhylp123/hifiasm)
- PacBio reads were mapped back to the genome assembly with minimap2
- Illumina HiC reads were mapped back to the genome assembly with bwa mem

- PacBio and HiC reads were extracted from the alignment bam files:
   - from a 2Mb region within a 28Mb contig
   - a 20Kb gap was introduces at Position 1Mb and all PacBio and HiC reads that mapped into that region were filtered out
   - HiC reads were further filtered down to 50%, representing a 38X coverage

