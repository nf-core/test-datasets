# TAPS test data

Test data for TAPS (TET-assisted pyridine borane sequencing) in nf-core/methylseq: `--aligner bwamem3 --taps` (`bwa-mem3 mem --meth=taps`) with rastair (`--taps_caller rastair` or `rastair2`), used by `-profile test_taps` and `tests/bwamem3_taps.nf.test`.

The bisulfite test reads on this branch can't be used for this: `SRR389222` is 36 bp, too short for bwa-mem3 to map, and the E. coli reads don't match `reference/genome.fa`.

## Files

| Path                                        | Contents                                                                                                                     |
| ------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------- |
| `testdata/taps/taps_mm10_chr19_R1.fastq.gz` | 12,417 TAPS read pairs, read 1 (also used alone as a single-end sample)                                                        |
| `testdata/taps/taps_mm10_chr19_R2.fastq.gz` | TAPS reads, read 2                                                                                                            |
| `reference/genome_taps.fa.gz`               | Reference: a 41,226 bp region of mouse (mm10) chr19 (contig `chr19`), the CpG-methylated lambda spike-in (`bacteriophage_lambda_CpG`) and an unmodified 2 kb control (`2kb_3_Unmodified`) |
| `reference/genome_taps.fa.fai`              | samtools faidx index of the uncompressed reference                                                                            |

## Source

The reads are a subset of the mouse E14 ES cell TAPS data from Liu et al. (2019), *Bisulfite-free direct detection of 5-methylcytosine and 5-hydroxymethylcytosine at base resolution*, Nature Biotechnology 37, 424–429, which is public in GEO as [GSE112520](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112520) (TAPS sample GSM3454335, SRA runs SRR8145389–SRR8145392). They were taken from the test data of [rastair](https://github.com/bsbludwig/rastair) v2.2.0 (`tests/data/`); its author, Benjamin Schuster-Böckler, confirmed on 2026-11-09 that the data can be used here.

- **Reads:** every primary read pair in rastair's `tests/data/test.bam`, converted to FASTQ:

  ```
  samtools view -u -F 0x900 test.bam | samtools collate -u -O - \
    | samtools fastq -1 taps_mm10_chr19_R1.fastq.gz -2 taps_mm10_chr19_R2.fastq.gz -s /dev/null -0 /dev/null -n -
  ```

- **Reference:** rastair's `tests/data/test.fasta.gz` with `chr19` cut down to mm10 chr19:6,102,575-6,143,800 (the region the reads come from, kept under the contig name `chr19`, so positions are relative to the start of that region), and the `bacteriophage_lambda_CpG` and `2kb_3_Unmodified` spike-ins unchanged. The full chr19 file is 19.5 MB; the cut-down one keeps the tests fast.

Expected results (bwa-mem3 `--meth=taps`, rastair 0.8.2 and 2.2.0): ~99.8% of reads mapped, CpG methylation of ~0.98 on the lambda spike-in, ~0.004 on the unmodified 2 kb control and ~0.38 on chr19.
