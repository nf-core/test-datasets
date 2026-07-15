# Synthetic trio data for MerquryFK

To test MerquryFK in trio mode, we need a true trio dataset to avoid floating-point errors during testing.

This dataset copies a (compressed) version of the synthetic trio dataset included in the MerquryFK repository: https://github.com/thegenemyers/MERQURY.FK/tree/main/SMALL_SYNTHETIC_TRIO

## Contents

- Child.Reads.fasta.gz - synthetic HiFi reads for the child in the trio
- Mother.Reads.fasta.gz - synthetic HiFi reads for the mother in the trio
- Father.Reads.fasta.gz - synthetic HiFi reads for the father in the trio
- Hap1.fasta.gz - assembled first haplotype of the child
- Hap2.fasta.gz - assembled second haplotype of the child
