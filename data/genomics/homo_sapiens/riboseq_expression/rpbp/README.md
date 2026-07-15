# Test data for `rpbp/*` modules

Per-stage intermediates from a single end-to-end `fasta_gtf_bam_rpbp` subworkflow run on the existing chr20 fixture (`SRX11780888_chr20.bam` + `Homo_sapiens.GRCh38.111_chr20.gtf`). Each fixture is the immediate-upstream input for one rpbp module, so module-level tests can fetch one static file rather than chain six upstream stages.

## Why these fixtures?

`modules/nf-core/rpbp/*` test setups used to chain `GUNZIP -> PREPAREGENOME -> EXTRACTMETAGENEPROFILES -> ESTIMATEMETAGENEBAYESFACTORS -> SELECTPERIODICOFFSETS -> GETPERIODICLENGTHSOFFSETS -> EXTRACTORFPROFILES -> ESTIMATEORFBAYESFACTORS` to test downstream modules like `SELECTFINALPREDICTIONSET`. Every chained run cost several minutes of CI time per module. With these fixtures, each module test fetches its one immediate-upstream output and runs in well under a minute.

The end-to-end integration test still lives in `subworkflows/nf-core/fasta_gtf_bam_rpbp/tests/main.nf.test`, which exercises the full chain.

## Files

| File                                                          | Size    | Description                                                                                                                                                                                              |
| ------------------------------------------------------------- | ------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `reference.annotated.bed.gz`                                  | <500 KB | Transcript-level annotated BED output by `rpbp/preparegenome`. Consumed by `rpbp/extractmetageneprofiles`.                                                                                               |
| `reference.orfs-genomic.annotated.bed.gz`                     | <500 KB | Genomic-coordinate ORF BED output by `rpbp/preparegenome`. Consumed by `rpbp/extractorfprofiles` and `rpbp/estimateorfbayesfactors`.                                                                     |
| `reference.orfs-exons.annotated.bed.gz`                       | <500 KB | Exon-coordinate ORF BED output by `rpbp/preparegenome`. Consumed by `rpbp/extractorfprofiles`.                                                                                                           |
| `SRX11780888_chr20.metagene-profile.csv.gz`                   | <50 KB  | Per-read-length metagene profile output by `rpbp/extractmetageneprofiles`. Consumed by `rpbp/estimatemetagenebayesfactors`.                                                                              |
| `SRX11780888_chr20.metagene-periodicity-bayes-factors.csv.gz` | <50 KB  | Per-read-length periodicity Bayes-factor table output by `rpbp/estimatemetagenebayesfactors`. Consumed by `rpbp/selectperiodicoffsets`.                                                                  |
| `SRX11780888_chr20.periodic-offsets.csv.gz`                   | <50 KB  | Per-read-length periodic-offset table output by `rpbp/selectperiodicoffsets`. Consumed by `rpbp/getperiodiclengthsoffsets`.                                                                              |
| `SRX11780888_chr20.periodic_lengths_offsets.tsv`              | <1 KB   | Filtered length/offset pairs output by `rpbp/getperiodiclengthsoffsets` using lenient `'10 1 None 0.0'` thresholds (chr20 alone does not pass the rpbp defaults). Consumed by `rpbp/extractorfprofiles`. |
| `SRX11780888_chr20.profiles.mtx.gz`                           | <2 MB   | Per-ORF Ribo-seq read-count profile sparse matrix output by `rpbp/extractorfprofiles`. Consumed by `rpbp/estimateorfbayesfactors`.                                                                       |
| `SRX11780888_chr20.bayes-factors.bed.gz`                      | <2 MB   | Per-ORF translation-Bayes-factor table output by `rpbp/estimateorfbayesfactors`. Consumed by `rpbp/selectfinalpredictionset`.                                                                            |

All files <4 MB. Total set <10 MB.

## How they were derived

A single `nf-core subworkflows test fasta_gtf_bam_rpbp` run, with input

- BAM: `aligned_reads/SRX11780888_chr20.bam` (+ `.bai`)
- FASTA: `Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz` (gunzipped)
- GTF: `Homo_sapiens.GRCh38.111_chr20.gtf`

from this same `riboseq_expression/` folder, plus the module-level test config:

```
process {
    withName: 'RPBP_GETPERIODICLENGTHSOFFSETS' {
        ext.args = '10 1 None 0.0'
    }
    withName: 'RPBP_SELECTFINALPREDICTIONSET' {
        ext.args = '--select-longest-by-stop --select-best-overlapping'
    }
}
```

The `RPBP_GETPERIODICLENGTHSOFFSETS` thresholds (`min_metagene_profile_count=10`, `min_metagene_bf_mean=1`, `max_metagene_bf_var=None`, `min_metagene_bf_likelihood=0.0`) are intentionally lenient so chr20 produces non-empty output; rpbp's production defaults are tuned for whole-genome data and reject chr20-only profiles.

Output files were captured directly from each module's work directory.

## rpbp version

`rpbp 4.0.1` (Wave container `community.wave.seqera.io/library/rpbp:4.0.1--71297b462026e13b`).

Used by `modules/nf-core/rpbp/*/tests/main.nf.test`.
