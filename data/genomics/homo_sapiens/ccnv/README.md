# copynumber-cnv-testdata

Test datasets and fixtures for [`ccnv`](https://github.com/matthdsm/copynumber-cnv).

## Data index

### Synthetic alignment fixtures
Minimal deterministic synthetic alignment fixtures (`136 KB` total) for end-to-end testing of `ccnv extract` across file formats, index types, contig naming conventions, read pairing modes, and filter flags.

| File | Format / Index | Layout / Build | Description |
| --- | --- | --- | --- |
| [`synthetic_ref.fa`](synthetic_ref.fa) | FASTA (`.fai`) | Synthetic UCSC (`chr1` 25 kb, `chr2` 12,345 bp, `chrX` 10 kb, `chrY` 10 kb, `chrM` 5 kb, `chrUn_decoy` 5 kb) | Minimal 67.3 kb reference genome with varying per-bin GC content, a partial trailing bin on `chr2`, sex contigs, and non-canonical contigs. |
| [`paired_end.bam`](paired_end.bam) | BAM (`.bai`) | Synthetic UCSC (6 contigs) | Paired-end proper-pair reads (`80` pairs total; `72` canonical, `8` non-canonical) spanning short (`[80, 150]` bp), mononucleosomal, dinucleosomal, and high-MW fragment bands; tests BAM + `.bai` extraction, multi-threaded tiling (`-t > 1`), fragmentomics, and `--all-contigs`. |
| [`paired_end.cram`](paired_end.cram) | CRAM (`.crai`) | Synthetic UCSC (6 contigs, ref `synthetic_ref.fa`) | Bit-for-bit identical reads to `paired_end.bam`; tests CRAM + `--fasta`, `.crai` tile streaming, and BAM/CRAM extraction parity. |
| [`paired_end.ccnv`](paired_end.ccnv) | `2.0.0` (sample) | Synthetic UCSC (6 contigs) | extracted from `paired_end.bam` |

### Synthetic reference cohort

Minimal deterministic cohort (`640 KB` total) for end-to-end testing of `ccnv reference`, `ccnv ref-qc`, samplesheet formats, rebinning, and schema rejection.

| File | Format / Schema | Layout / Mode | Key Properties | Description |
| --- | --- | --- | --- | --- |
| [`samples/`](samples/) | `2.0.0` (12 samples) | Synthetic 6-contig | 5 kb bins, paired-end, 6 XX (`ctrl_f01`–`f06`) + 6 XY (`ctrl_m01`–`m06`) | Minimal balanced control cohort for `ccnv reference` end-to-end tests, 2D sex inference, and covariate normalization. |
| [`samplesheet_germline.yaml`](samplesheet_germline.yaml) | YAML | 12 entries (unlabeled `sex`) | Relative paths to `samples/` | Tests YAML samplesheet loading and automated 2D sex inference. |
| [`ref_qc_warn_flat.ccnv`](ref_qc_warn_flat.ccnv) | `2.0.0` (reference) | Synthetic UCSC (`synthetic_ref.fa`) / flat | 5 kb bins, `0` samples (`flat=true`) | Flat FASTA reference fixture verified to produce `ref-qc` verdict `WARN`. |