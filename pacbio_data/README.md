# PacBio Test Data

This directory contains PacBio test datasets used for local pipeline testing.

## `revio-with-m6a-tags.custom-reference.intervals.bed`

This BED file remaps the five regions in `revio-with-m6a-tags.intervals.bed` to
the contigs in the reduced reference described below. Each reduced contig has
25,000 bp of reference sequence on both sides of the original interval, so the
target begins at zero-based offset 25,000. The file was generated with:

```bash
awk 'BEGIN { OFS="\t" }
{
    reference_start = $2 - 24999
    reference_end = $3 + 25000
    print $1 "_" reference_start "_" reference_end, 25000, 25000 + ($3 - $2)
}' revio-with-m6a-tags.intervals.bed \
    > revio-with-m6a-tags.custom-reference.intervals.bed
```

The subtraction uses 24,999 rather than 25,000 because BED starts are
zero-based, whereas the region starts passed to `samtools faidx` are one-based.

## `GRCh38_revio_m6a_intervals.fasta` and `GRCh38_revio_m6a_intervals.fasta.fai`

This 375,616-bp reduced reference was curated from the GATK GRCh38 iGenomes
reference, `Homo_sapiens_assembly38.fasta`. It contains the five regions in
`revio-with-m6a-tags.intervals.bed`, expanded by 25,000 bp on each side. This
avoids staging the full GATK GRCh38 reference and indexes during lightweight
PacVar tests, reducing both runtime and disk usage.

The custom FASTA was made by adding 25,000 bp of flanking sequence to each
interval, extracting those five regions from the GATK GRCh38 reference, naming
each resulting contig for its original chromosome and coordinates, and indexing
the reduced reference with `samtools faidx`.

## `revio-with-kinetics.intervals.bed`
An compact BED file covering the long-read alignments in the Revio Fiberseq kinetics test BAM. These intervals reduce DeepVariant runtime in the PacVar Fiber-seq kinetics test profile that uses `revio-with-kinetics.bam` as the test set and `test_fiberseqs_with-kinetics_tags.config` as the test config.

## `revio-with-with-m6A-tags.intervals.bed`
An compact BED file covering the long-read alignments in the Revio Fiberseq m6A test BAM. These intervals reduce DeepVariant runtime in the PacVar Fiber-seq with m6A tags test profile that uses `revio-with-m6a-tags.bam` as the test set and `test_fiberseqs_with_m6A_tags.config` as the test config.



## `revio-with-kinetics.bam`

`revio-with-kinetics.bam` is a PacBio long-read Fiber-seq test dataset, downsized (10 reads)
from the[`fiberseq/fibertools-rs`](https://github.com/fiberseq/fibertools-rs)
repository. It contains kinetic signatures and no m6A tags (A+a). 


Upstream source:
[`tests/data/rebio.bam`](https://github.com/fiberseq/fibertools-rs/blob/main/tests/data/rebio.bam)

## `revio-with-m6a-tags.bam`

`revio-with-m6a-tags.bam` is a PacBio long-read Fiber-seq 100-read subset of the
source dataset in GSE330647, obtaining m6A modification tags (A+a).

