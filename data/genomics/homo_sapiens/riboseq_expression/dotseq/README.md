# Test data for `dotseq/dotseq`

A small cohort of ORF-level Ribo-seq + RNA-seq counts with matched ORF annotation, derived from the `cell_cycle_subset` example data shipped in the [DOTSeq Bioconductor package](https://bioconductor.org/packages/release/bioc/html/DOTSeq.html) (`inst/extdata`).

## Why this set?

The `dotseq/dotseq` module wraps `DOTSeqDataSetsFromSummarizeOverlaps()`, which takes a per-ORF count matrix plus a `GRanges` ORF annotation. These fixtures provide both in tidy TSV form (one ORF per row), aligned on a stable `orf_id` so the module's R template can rebuild the `GRanges` and run end-to-end against a known-good cohort the package author already validated against the API.

## Files

| File | Size | Description |
|---|---|---|
| `counts.tsv.gz` | 73 KB | Per-ORF count matrix (6642 ORFs x 12 samples). First column is `orf_id`, remaining columns are sample IDs (6 Ribo-seq + 6 RNA-seq). Sample IDs match the `run` column of `samplesheet.csv`. |
| `annotation.tsv.gz` | 72 KB | Per-ORF annotation (6642 rows). Columns: `orf_id, gene_id, chrom, start, end, strand, orf_type`. `orf_type` is one of `mORF`, `uORF`, `dORF`. |
| `samplesheet.csv` | 423 B | 12-sample condition table covering 3 conditions (Mitotic_Cycling, Mitotic_Arrest, Interphase) x 2 strategies (ribo, rna). Columns: `run, strategy, replicate, condition`. |
| `metadata.txt.gz` | 211 B | DOTSeq's headerless 24-sample metadata covering both chx + har treatment arms, columns: `run strategy replicate treatment condition`. Kept verbatim from the upstream package for traceability; not consumed by the nf-test directly. |

## How they were derived

The fixtures are derived from `system.file("extdata", package = "DOTSeq")` in the Bioconductor 3.23 release of DOTSeq (v1.0.0), itself sourced from the `compgenom/DOTSeq` GitHub repository at the same tag. Upstream cohort: GSE231096 cell-cycle Ribo-seq + RNA-seq (Ly et al. 2024), restricted to the chx (cycloheximide) treatment arm.

- `counts.tsv.gz`: extracted from `featureCounts.cell_cycle_subset.txt.gz` by dropping the featureCounts header banner and the 5 annotation columns (`Chr, Start, End, Strand, Length`), renaming `Geneid` to `orf_id` (with a within-gene `:O<###>` suffix that matches DOTSeq's internal ORF naming), and stripping the BAM-path sample columns down to their SRR accessions.
- `annotation.tsv.gz`: joined from `gencode.v47.orf_flattened_subset.gtf.gz` (genomic coordinates + `gene_id`) and `gencode.v47.orf_flattened_subset.bed.gz` (`orf_type` per ORF), keyed on `chrom + start + end + strand + gene_id`. ORF ids match `counts.tsv.gz`.
- `samplesheet.csv`: built from `metadata.txt.gz` by adding column headers (`run, strategy, replicate, treatment, condition`) and filtering to `treatment == "chx"` (matching the 12 sample columns the count table actually contains). The redundant `treatment` column is then dropped.

## Source & licence

[compgenom/DOTSeq](https://github.com/compgenom/DOTSeq), MIT licence. Upstream sample accessions: `SRR24230462` to `SRR24230485` (GSE231096, Ly et al. 2024).

Used by `modules/nf-core/dotseq/dotseq/tests/main.nf.test`.
