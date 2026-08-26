# nf-core/dartseq Test Dataset

## Overview

This directory contains test data for the nf-core/dartseq pipeline, which detects and quantifies m6A RNA modifications from DART-seq data by identifying APOBEC1-YTH-induced C-to-U editing events adjacent to modified adenosines.

Two kinds of test data are provided:

- A small, fully synthetic 4-sample dataset (`fastq/`, `samplesheet_4samples.csv`, `reference/`) used by the pipeline's default `-profile test`, covering the full alignment + Bullseye site-calling workflow end-to-end.
- A subsetted real-data fixture (`glm/`) derived from an actual DART-seq experiment, used specifically to exercise the `BULLSEYE_R_GLM` differential-editing statistics step with real coverage/mutation counts.

## Files

### Synthetic 4-sample dataset

- `fastq/{A03,A04,C03,E01}_R{1,2}.fastq.gz` - Paired-end DART-seq reads for 4 samples (25,000 read pairs each). A03/C03 are `control` group, A04/E01 are `dart` (edited) group.
- `samplesheet_4samples.csv` - nf-core-format samplesheet (`sample,fastq_1,fastq_2,group`) referencing the fastqs above by their `raw.githubusercontent.com` URL.
- `contrasts_test.csv` - Example custom contrasts file (`contrast_id,edited_group,control_group,mode,min_edit,max_edit,fold_threshold,min_sites`) for testing the `--bullseye_contrasts` option.
- `reference/genome.fa` - Single-contig (`chr22_test`) reference genome, 180,001 bp, covering an ~11-gene region of human chr22 (DDT, DDTL, DERL3, GSTT2, GSTT2B, GSTTP2, LOC391322, MIF, MIF-AS1, SLC2A11, SMARCB1).
- `reference/genes.gtf` - GTF annotation matching `genome.fa` (286 records, same gene set as above).
- `reference/bullseye.test.refFlat` - refFlat-format annotation matching `genome.fa`/`genes.gtf`, used by Bullseye's `Find_edit_sites.pl`/`quantify_sites.pl` (23 transcripts).

### GLM real-data fixture (`glm/`)

- `glm/coverage.txt`, `glm/mut.txt`, `glm/score.txt` - `BULLSEYE_GATHER_SITES` output matrices (per-site coverage / mutation count / score), restricted to `chr19` (3,143 sites). Derived from a single real DART-seq A04-vs-A03 contrast, reprocessed through the current pipeline's `BULLSEYE_QUANTIFY_SITES` + `BULLSEYE_ACFILTER` + `BULLSEYE_GATHER_SITES` modules to guarantee the file format matches what the pipeline currently produces.
- `glm/observations.tsv` - Matching `BULLSEYE_R_GLM` design/colData file (`sample`, `group` columns).

  Only one real sample/contrast is available in this fixture, so the beta-binomial GLM design has no residual degrees of freedom - `glm_comp_bbn()` in Bullseye's `functions.R` returns all-`NA` p-values by design whenever a comparison factor has fewer than 2 levels, rather than erroring. This fixture exercises `BULLSEYE_R_GLM`'s file-format handling and R execution path with real biological data; it is not a statistically meaningful GLM comparison (that would require multiple real contrasts, which this dataset does not have).

**Data provenance and permission:** the `glm/` files derive from RNA-seq data generated for a METTL3 DART-seq study; redistributed here in subsetted, aggregated form (per-site coverage/mutation counts only - no raw reads or alignments) with permission from the data owner.

## Usage

Use the test profile in nf-core/dartseq:

```bash
nextflow run nf-core/dartseq -profile test,docker
```

The test configuration automatically references the samplesheet and reference files from this repository via `params.pipelines_testdata_base_path`.
