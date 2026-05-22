# Test data for `dotseq/dotseq`

The bundled `cell_cycle_subset` example data from the [DOTSeq Bioconductor package](https://bioconductor.org/packages/release/bioc/html/DOTSeq.html), copied verbatim from its `inst/extdata` plus a derived headered samplesheet. Total: ~252 KB across 5 files.

## Why this set?

DOTSeq's `DOTSeqDataSetsFromFeatureCounts()` requires four inputs that share an internal contract: an ORF-level featureCounts table (`Geneid, Chr, Start, End, Strand, Length, samples...`), a flattened ORF GTF (with `gene_id` + `exon_number` attributes), a matching flattened BED, and a condition table with `run, strategy, replicate, condition` columns where the ORF naming and sample IDs all line up. The existing `riboseq_expression/` fixtures provide gene-level Salmon counts + a regular Ensembl GTF; neither matches DOTSeq's contract. Re-deriving compatible inputs would require running an ORF caller + featureCounts on the existing BAMs. Instead this PR ships the small cohort the package author already validated against the API, which makes the module test reproducible against a known-good fixture set.

## Files

| File | Size | Description |
|---|---|---|
| `featureCounts.cell_cycle_subset.txt.gz` | 117 KB | featureCounts v2.1.1 output, 6644 ORF rows × 12 sample columns. Subset of GSE231096 cell-cycle Ribo-seq + RNA-seq filtered to the chx (cycloheximide) treatment arm: 6 Ribo-seq + 6 RNA-seq samples across 3 conditions (Mitotic_Cycling, Mitotic_Arrest, Interphase). |
| `gencode.v47.orf_flattened_subset.gtf.gz` | 81 KB | Flattened GENCODE v47 ORF annotation (6945 lines), matching the count table on `gene_id:O###` naming. |
| `gencode.v47.orf_flattened_subset.bed.gz` | 53 KB | Same ORFs as the GTF in BED format (6642 lines). |
| `metadata.txt.gz` | 211 B | DOTSeq's headerless 24-sample metadata covering both chx + har (harringtonine) treatments, columns: `run strategy replicate treatment condition`. Kept verbatim from the package. |
| `samplesheet.csv` | 423 B | Headered, chx-only subset of `metadata.txt.gz` (12 rows: 6 Ribo + 6 RNA) - what the nf-test consumes directly. |

## How they were derived

1. The first four files come straight from `system.file("extdata", package = "DOTSeq")` in the Bioconductor 3.23 release of DOTSeq (v1.0.0), which is itself sourced from the `compgenom/DOTSeq` GitHub repository at the same tag.
2. `samplesheet.csv` is built from `metadata.txt.gz`: column headers `run,strategy,replicate,treatment,condition` were added, and rows were filtered to `treatment == "chx"` (matching the 12 sample columns the bundled featureCounts table actually contains). The `treatment` column is then redundant and dropped to keep the four-column shape DOTSeq's `parse_condition_table()` expects.

## Verified

Running `nf-core modules test --profile docker dotseq/dotseq` against this fixture set with contrast Mitotic_Cycling vs Interphase produces a complete output set in ~4 min on a c5.9xlarge: per-ORF DTE + DOU contrast tables (`translation.dotseq.results.tsv`, `dou.dotseq.results.tsv`), per-condition strategy contrasts, the four `plotDOT()` PNGs (volcano / composite / venn / heatmap), a DTE p-value distribution histogram, the serialised `DOTSeqDataSets.rds`, R sessionInfo and versions.yml.

## Source & licence

[compgenom/DOTSeq](https://github.com/compgenom/DOTSeq), MIT licence. Upstream sample accessions: `SRR24230462`–`SRR24230485` (GSE231096, Ly et al. 2024).

Used by `modules/nf-core/dotseq/dotseq/tests/main.nf.test`.
