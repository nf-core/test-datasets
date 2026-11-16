# nf-core/abotyper test fixtures

Small files used by the nf-test tests of nf-core/abotyper, so each step can be tested without running
the steps before it.

Per-sample intermediate files are in `../per_sample_processing/<sample>/`, laid out like the pipeline's
own output (`alignment/`, `clair3/`, `clair2metrics/`, `mpileup/` and the sample's `.ABOPhenotype.txt`).
They come from the pipeline's test run on `ABO_AA_barcode01` and `ABO_BO1_barcode06` (the FASTQ files
in `../`), aligned to the full-length ABO RefSeqGene NG_006669.2.

This folder holds the rest:

- `samples_two.tsv`, `samples_named.tsv`, `controls_AA_negative.tsv`: sample tables for
  `abotyper-aggregate --samples`.
- `empty.fastq.gz`: a FASTQ with no reads, used as a clean negative control.
- `samplesheet_controls_multirow.csv`: samplesheet for the second pipeline test (a zero-read negative
  control, a numeric barcode and a sample with three rows whose files share a name).
