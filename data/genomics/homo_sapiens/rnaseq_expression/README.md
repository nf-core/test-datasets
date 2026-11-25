# rnaseq_expression

Exon and transcript count tables for four samples of human chromosome X from the Geuvadis project: two GBR (`ERR188383`, `ERR188428`) and two YRI (`ERR188454`, `ERR204916`). All inputs come from the `rnasplice` branch. Each table is cut down to a few genes with more than one feature to keep it small.

- `chrX.samplesheet.csv`: sample IDs (`sample`) and population (`condition`).
- `chrX.exon_counts.tsv`: exon counts with `gene_id`, `exon_id` and one column per sample. Made with `dexseq/count` (`--order pos --stranded no --minaqual 10`) from `testdata/<sample>.Aligned.out.bam` and `testdata/dexseq/genes_chrX.DEXSeq.gff.gz`, merged into one table, and limited to the 14 most-expressed genes with 2 to 25 counted exonic parts.
- `chrX.exon_counts.featureCounts.tsv`: exonic parts counted with featureCounts v2.1.1 (`-F GTF -t exonic_part -g gene_id -f -O --countReadPairs -p -s 0`) from `chrX_subset/testdata/<sample>.Aligned.out.bam` and the same flattened annotation, in featureCounts' multi-sample format: a `#` command line, then `Geneid`, `Chr`, `Start`, `End`, `Strand`, `Length` and one column per BAM file. Limited to the 7 genes with counts.
- `chrX.transcript_counts.tsv`: Salmon estimated counts from `testdata/salmon_quant/salmon_quant.tar.gz`, imported with `tximeta/tximport` using a transcript-to-gene map built by `custom/tx2gene` from `reference/genes_chrX.gtf`. Columns are transcript ID (`tx`), gene ID (`gene_id`) and one column per sample, limited to the 40 most-expressed genes with at least two expressed transcripts.

Used by `edger/diffsplice`.
