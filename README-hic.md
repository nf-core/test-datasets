# Test datatest for nf-core-hic pipeline

Data From GSE87311 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87311)

See Schalbetter SA, Goloborodko A, Fudenberg G, Belton JM et al. SMC complexes differentially compact mitotic chromosomes according to genomic context. Nat Cell Biol 2017 Sep;19(9):1071-1080.
https://www.ncbi.nlm.nih.gov/pubmed/28825700

This small dataset contains 500 000 reads from the G1-R1 samples
* Sample GSM2327656 (G1-R1)
* Saccharomyces cerevisiae
* strain: W303-1a
* Restriction enzyme: HindIII

Command line:
```
nextflow run main.nf --fasta ./reference/W303_SGD_2015_JRIU00000000.fsa --reads './data/*_R{1,2}.fastq.gz' --skip_cool
```