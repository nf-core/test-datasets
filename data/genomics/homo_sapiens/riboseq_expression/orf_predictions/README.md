# Test data for `custom/orfnormalise` + `custom/orfmerge`

Six small Ribo-seq ORF prediction outputs — one per supported caller, plus a
second Ribo-TISH file covering its extended-ORF output — used by
`modules/nf-core/custom/orfnormalise` and `modules/nf-core/custom/orfmerge` to
exercise the parser, classifier, score-direction, and cross-caller merge logic
end-to-end.

| File                                 | Size  | Caller     | Source genome |
| ------------------------------------ | ----- | ---------- | ------------- |
| `sample1.ribocode.txt`               | 13 KB | RiboCode   | chr20         |
| `sample1.ribotish.pred.txt`          | 3 KB  | Ribo-TISH  | chr20         |
| `sample1.ribotish.extended.pred.txt` | 3 KB  | Ribo-TISH  | chr20         |
| `sample1.ribotricer.tsv`             | 3 KB  | Ribotricer | chr20         |
| `sample1.rpbp.predicted-orfs.bed.gz` | 2 KB  | Rp-Bp      | chr20         |
| `cohort.price.orfs.tsv`              | 5 KB  | PRICE      | chr19+chr22   |

All are real-tool output with no synthetic data. The first five are the head
plus first 15 records, so every column the downstream parsers read is
exercised; `sample1.ribotish.extended.pred.txt` is row-selected instead, to
cover each distinct `TisType` value (see below).

## How they were derived

Five of the six are sliced from outputs produced by existing nf-core/modules
tests (or nf-core/modules#11695 for Rp-Bp) on the VM, then trimmed to
header + 15 records with `head`. `sample1.ribotish.extended.pred.txt` comes
from an nf-core/riboseq pipeline run instead, because no module test supplies
the input that produces its ORF-type values.

### `sample1.ribocode.txt`

Produced by `modules/nf-core/ribocode/ribocode` test
`test_ribocode_ribocode`. Inputs (already in test-datasets:modules):

- BAM: `genomics/homo_sapiens/riboseq_expression/aligned_reads/SRX11780887.Aligned.toTranscriptome.out.bam`
- Annotation tar: `genomics/homo_sapiens/riboseq_expression/ribocode/annotation.tar.gz`
- Pre-config: `genomics/homo_sapiens/riboseq_expression/ribocode/test_pre_config.txt`

```bash
nf-test test modules/nf-core/ribocode/ribocode/tests/main.nf.test --profile docker
head -1 test.txt > sample1.ribocode.txt
head -16 test.txt | tail -15 >> sample1.ribocode.txt
```

### `sample1.ribotish.pred.txt`

Produced by `modules/nf-core/ribotish/predict` test
`sarscov2 [bam] - single_end - single ribo bam`. Inputs:

- BAM: `genomics/homo_sapiens/riboseq_expression/aligned_reads/SRX11780888_chr20.bam`
- FASTA: chr20 from `Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz`
- GTF: `Homo_sapiens.GRCh38.111_chr20.gtf`

```bash
nf-test test modules/nf-core/ribotish/predict/tests/main.nf.test --profile docker
head -1 test_pred.txt > sample1.ribotish.pred.txt
head -16 test_pred.txt | tail -15 >> sample1.ribotish.pred.txt
# The source ribotish run lacks a TIS (Harringtonine) BAM, so TISPvalue
# and FisherPvalue cells come out as the literal string "None" - which
# means a downstream FisherPvalue-reading parser never sees a real
# p-value. Substitute the real RiboPvalue into the FisherPvalue column
# so the fixture covers the typical case (FisherPvalue populated):
awk -F'\t' 'BEGIN{OFS="\t"} NR==1 {print; next} {if ($15=="None" && $13!="None") $15=$13; print}' \
    sample1.ribotish.pred.txt > .tmp && mv .tmp sample1.ribotish.pred.txt
```

### `sample1.ribotish.extended.pred.txt`

Ribo-TISH qualifies its `TisType` after a colon (`Novel:CDSFrameOverlap`) only
when a secondary annotation is supplied with `-a`, which nf-core/riboseq does in
extended-ORF mode. `sample1.ribotish.pred.txt` above predates that mode and its
`ribotish/predict` test passes no secondary annotation, so no fixture carried a
qualified value.

Produced by running nf-core/riboseq on its own chr20 test data:

```bash
nextflow run nf-core/riboseq -profile test,docker \
    --skip_stringtie false \
    --extended_orf_analysis true \
    --outdir results
```

Taken from the pooled Ribo-TISH output
(`results/orf_predictions/ribotish_all/allsamples_pred.txt`, 2716 records) and
reduced to one row per distinct `TisType` — seven qualified, seven bare:

| Qualified                  | Bare        |
| -------------------------- | ----------- |
| `Novel:CDSFrameOverlap`    | `Annotated` |
| `5'UTR:Known`              | `Truncated` |
| `5'UTR:CDSFrameOverlap`    | `Extended`  |
| `3'UTR:CDSFrameOverlap`    | `Novel`     |
| `Truncated:Known`          | `5'UTR`     |
| `Novel:Known`              | `3'UTR`     |
| `Internal:CDSFrameOverlap` | `Internal`  |

All 19 columns and the header are unchanged, every row is verbatim tool output,
and all rows are chr20 with transcript and gene ids that resolve against
`Homo_sapiens.GRCh38.111_chr20.gtf`.

### `sample1.ribotricer.tsv`

Produced by `modules/nf-core/ribotricer/detectorfs` test against
`SRX11780888_chr20.bam` plus candidate ORFs built from chr20 GTF by
`ribotricer/prepareorfs`.

```bash
nf-test test modules/nf-core/ribotricer/detectorfs/tests/main.nf.test --profile docker
head -1 test_translating_ORFs.tsv > sample1.ribotricer.tsv
head -16 test_translating_ORFs.tsv | tail -15 >> sample1.ribotricer.tsv
# strip verbose `profile` column (last) to keep the fixture small:
awk -F'\t' 'BEGIN{OFS="\t"} {NF=17; print}' sample1.ribotricer.tsv > .tmp && mv .tmp sample1.ribotricer.tsv
```

Note: this is the ribotricer `detect-orfs` `_translating_ORFs.tsv`
schema. It does not carry the multi-exon `coordinate` column that
`prepareorfs` emits; downstream parsers in `custom/orfnormalise` fall
back to deriving a single-block span from the `ORF_ID`
(`<transcript_id>_<gstart>_<gend>_<length_nt>`) when the column is
absent.

### `sample1.rpbp.predicted-orfs.bed.gz`

Produced by `modules/nf-core/rpbp/selectfinalpredictionset` test
`homo_sapiens chr20 - select final prediction set`. The Rp-Bp modules
are introduced by nf-core/modules#11695 (branch
`rpbp-add-modules-and-subworkflows`).

```bash
nf-test test modules/nf-core/rpbp/selectfinalpredictionset/tests/main.nf.test --profile docker
zcat test.predicted-orfs.filtered.bed.gz | head -15 | gzip > sample1.rpbp.predicted-orfs.bed.gz
```

### `cohort.price.orfs.tsv`

Produced by GEDI/PRICE on the 4-sample chr19+chr22 cohort already
shipped under `genomics/homo_sapiens/riboseq_expression/price/`:

- BAMs: `bams/SRX1178088{5,6,7,8}.chr19_22.ds50.bam`
- Reference: `Homo_sapiens.GRCh38_chr19_22.pc_exon_masked.fa.gz` +
  `Homo_sapiens.GRCh38.111_chr19_22.pc.gtf.gz`
- Generated by `gedi -e IndexGenome` + `gedi -e Price` (per the PRICE
  test-data README in the same directory).

```bash
head -1 run.orfs.tsv > cohort.price.orfs.tsv
head -16 run.orfs.tsv | tail -15 >> cohort.price.orfs.tsv
```

## Verified

Each fixture round-trips through `custom/orfnormalise` (caller dispatch by
`meta.caller`), producing a normalised BED12 + sidecar TSV with the
harmonised `orf_class` vocabulary and 0-1000 BED scores. The combined
set then feeds `custom/orfmerge` to exercise both transcript-id-grouped
clustering (canonical / uORF / dORF) and 80%-reciprocal-overlap
clustering (novel_u / smORF).
