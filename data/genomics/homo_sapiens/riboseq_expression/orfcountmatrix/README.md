# Test data for `custom/orfcountmatrix`

Hand-crafted synthetic fixtures for the `custom/orfcountmatrix` module, which pivots per-sample ORF P-site count TSVs into a single ORF x sample matrix. The data is small enough to be read at a glance and is designed so the expected output exercises every code path:

- ORFs present in both samples (counts pulled from both columns).
- ORFs present in only one of the two samples (zero-fill in the other).
- ORFs present in the catalogue but in neither sample (row of zeros).

No real Ribo-seq data is needed here - the module is a deterministic pivot of integers, not a biology tool.

## Files

| File | Size | Description |
|---|---|---|
| `orf_catalogue.bed12` | 1.1 KB | 15-ORF BED12 catalogue. Column 4 is the ORF id, which defines row order and membership of the output matrix. |
| `samples/sample1.orf_psite_counts.tsv` | 249 B | 8-row count TSV for sample1 (`sample_id`, `orf_id`, `count`). Covers catalogue ORFs 1-8. |
| `samples/sample2.orf_psite_counts.tsv` | 250 B | 8-row count TSV for sample2 (`sample_id`, `orf_id`, `count`). Covers catalogue ORFs 1-3 and 9-13. |

Total: ~1.6 KB across 3 files.

## Coverage

- Catalogue ORFs 1-3: in both samples (intersection).
- Catalogue ORFs 4-8: sample1 only.
- Catalogue ORFs 9-13: sample2 only.
- Catalogue ORFs 14-15: catalogue-only (absent from both samples; should appear as a row of zeros in the output matrix).

## Expected output

`orf_psite_counts.tsv` with header `orf_id\tsample1\tsample2`, 15 data rows in catalogue order, with zeros for samples missing an ORF.

Used by `modules/nf-core/custom/orfcountmatrix/tests/main.nf.test`.
