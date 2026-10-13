# dexseq/exon

Exon counts for four samples of human chromosome X: two GBR (`ERR188383`, `ERR188428`) and two YRI (`ERR188454`, `ERR204916`), from the Geuvadis project. The counts were made with `dexseq/count` (`--order pos --stranded no --minaqual 10`) from `testdata/<sample>.Aligned.out.bam` and `testdata/dexseq/genes_chrX.DEXSeq.gff.gz` of the `rnasplice` branch.

- `<sample>.txt`: exon counts per sample, in the format written by `dexseq/count`. Only exonic parts with a non-zero count in at least one sample are kept, and only genes with at least two of these parts.
- `chrX.DEXSeq.gff.gz`: the flattened annotation, restricted to the same genes and exonic parts.
- `chrX.samplesheet.csv`: sample IDs (`sample`) and population (`condition`).

Used by `dexseq/exon`.
