# TAPS test data

Test data for TAPS (TET-assisted pyridine borane sequencing) in nf-core/methylseq: `--aligner bwamem3` (`bwa-mem3 mem --meth=taps`) with rastair (`--taps_caller rastair` or `rastair2`), used by `-profile test_taps` and `tests/bwamem3_taps.nf.test`.

The bisulfite data elsewhere on this branch cannot stand in for TAPS: in TAPS mode, bwa-mem3 maps almost none of those reads (bisulfite converts unmethylated rather than methylated cytosines).

## Files

| Path                                       | Contents                                                                                          |
| ------------------------------------------ | ------------------------------------------------------------------------------------------------- |
| `testdata/taps/taps_mm10_chr19_R1.fastq.gz` | TAPS reads, read 1 (also used alone as a single-end sample)                                        |
| `testdata/taps/taps_mm10_chr19_R2.fastq.gz` | TAPS reads, read 2                                                                                 |
| `reference/genome_taps.fa.gz`              | Reference: a ~41 kb region of mouse chr19 (contig `chr19`), plus the CpG-methylated lambda and unmodified 2 kb spike-ins |
| `reference/genome_taps.fa.fai`             | samtools faidx index of the uncompressed reference                                                 |

Every file must be redistributable under this repository's MIT licence. Each file should stay well under 5 MB.

Expected results (bwa-mem3 `--meth=taps`, rastair 0.8.2 and 2.2.0): ~99.8% of reads mapped, CpG methylation of ~0.98 on lambda, ~0.004 on the unmodified 2 kb control and ~0.38 on chr19.
