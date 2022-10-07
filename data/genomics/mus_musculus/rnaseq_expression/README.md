# Test Data

# SRP254919.salmon.merged.gene_counts.top1000cov.tsv

This data is intended for use in expression analysis downstream of the nf-core RNA-seq workflow. The objective is to make subset sufficiently to make the dataset manageable, but select genes in such a way that the data is still useful (e.g. we don't have a bunch of all-0 rows). [SRP254919](https://www.ebi.ac.uk/ena/browser/view/PRJNA622544?show=reads) is a small 6-sample mouse expression dataset that seemed a good candidate, and the data has been processed as follows: 

 * Sample sheet and FASTQ files derived by running [fetchngs](https://nf-co.re/fetchngs) with accession SRP254919.
 * Count matrix derived by running [rnaseq](https://nf-co.re/rnaseq) against the mouse genome GRCm38 with default parameters.
 * The top 1000 most variable genes from `salmon.merged.gene_counts` as follows: 

```
expression <- read.table("salmon.merged.gene_counts.tsv")

# Remove rows with zeros
expression <- expression[apply(expression[,c(-1,-2)], 1, function(x) all(x > 0)), ]

# Sort by coefficient of variation
expression <- expression[order(apply(expression[,c(-1,-2)], 1, function(x) sd(x)/mean(x)), decreasing = TRUE), ]

write.table(expression[1:1000,], file = "SRP254919.salmon.merged.gene_counts.top1000cov.tsv", sep="\t", row.names = FALSE, quote = FALSE)
```
