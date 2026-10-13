# Human RNA-seq expression

## chrX transcript counts

Transcript counts for four samples of human chromosome X: two GBR (`ERR188383`, `ERR188428`) and two YRI (`ERR188454`, `ERR204916`), from the Geuvadis project. The counts were imported with `tximeta/tximport` from the Salmon results in `testdata/salmon_quant/salmon_quant.tar.gz` of the `rnasplice` branch, with a transcript-to-gene map built by `custom/tx2gene` from `reference/genes_chrX.gtf` of the same branch.

- `chrX.transcript_counts.tsv`: transcript ID (`tx`), gene ID (`gene_id`) and estimated counts per sample. Only transcripts with a non-zero count in at least one sample are kept, and only genes with at least two of these transcripts, so that every gene can be tested for differential transcript usage.
- `chrX.samplesheet.csv`: sample IDs (`sample`) and population (`condition`).

Used by `dexseq/dtu`.
