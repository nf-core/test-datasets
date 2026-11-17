# Plasmodium falciparum WGS test reads

The six paired-end Illumina FASTQ files listed below are small, single-interval
subsets of whole-genome sequencing (WGS) data, not whole-genome coverage
datasets or amplicon reads.

## Samples and filenames

The six files in the table were copied without recompression from `ci_fastq/` in
the source dataset and constitute the documented CI dataset.

| Run accession | Source files in `ci_fastq/`  | Repository files in this directory        | Read pairs |
| ------------- | ---------------------------- | ----------------------------------------- | ---------: |
| ERR10375861   | `ERR10375861_{1,2}.fastq.gz` | `ERR10375861_Pf3D7_01_v3_R{1,2}.fastq.gz` |      9,579 |
| ERR10375871   | `ERR10375871_{1,2}.fastq.gz` | `ERR10375871_Pf3D7_01_v3_R{1,2}.fastq.gz` |        770 |
| ERR10375886   | `ERR10375886_{1,2}.fastq.gz` | `ERR10375886_Pf3D7_01_v3_R{1,2}.fastq.gz` |      1,876 |

In each row, source mate `_1` maps to repository `_R1`, and `_2` to `_R2`.
The source samplesheet identifies all three runs as `WGS`, `ILLUMINA`.

## Reference and interval

- Reference: [`../reference/Pf3D7_chr1.fa`](../reference/Pf3D7_chr1.fa),
  containing the complete `Pf3D7_01_v3` chromosome (640,851 bp).
- BED: [`../reference/Pf3D7_chr1.bed`](../reference/Pf3D7_chr1.bed),
  containing the single interval `Pf3D7_01_v3 184999 205000`.
- BED coordinates are **0-based, half-open**: bases 185000–205000 inclusive
  in 1-based coordinates (20,001 bp).

The extraction fetches mates, so reads outside the selected interval can be
included. The reference is the full chromosome, not a cropped interval.

## Provenance and extraction

The immediate source was the local dataset
`/home/vmduarut/Projects/pathogenepidemiology/Pf3D7_01_v3/`:
its `ci_fastq/` files, `Pf3D7_01_v3_CI.bed`, `samplesheet_Pf3D7_01_v3.csv`,
`README.md`, and `downsample.sh`. The chromosome reference at
`/home/vmduarut/Projects/pathogenepidemiology/Pf3D7_chr1.fa` is byte-identical
to the existing repository reference, which was left unchanged. No reference
indexes are stored alongside the repository FASTA.

The source README describes extraction from three coordinate-sorted, indexed
whole-genome BAMs aligned to an indexed whole-genome FASTA. The source extraction
script uses the same mate-fetching and FASTQ conversion workflow. The following
adaptation writes the repository-compatible names into a **new** output
directory; it requires Bash, samtools with `view -P`, and gzip. Set the reference
and BAM paths to the original inputs, and ensure the output directory's parent
exists. A failed run may leave partial output; inspect it before retrying.

```bash
set -euo pipefail
reference='/path/to/whole-genome.fa'
BAM_ERR10375861='/path/to/ERR10375861.bam'
BAM_ERR10375871='/path/to/ERR10375871.bam'
BAM_ERR10375886='/path/to/ERR10375886.bam'
out='/path/to/new-chr1-dataset'

mkdir -- "$out"
out=$(cd -- "$out" && pwd -P)
printf 'Pf3D7_01_v3\t184999\t205000\n' > "$out/Pf3D7_chr1.bed"
samtools faidx "$reference" Pf3D7_01_v3 > "$out/Pf3D7_chr1.fa"
samtools faidx "$out/Pf3D7_chr1.fa"

declare -A bams=(
  [ERR10375861]="$BAM_ERR10375861"
  [ERR10375871]="$BAM_ERR10375871"
  [ERR10375886]="$BAM_ERR10375886"
)
for sample in ERR10375861 ERR10375871 ERR10375886; do
  prefix="${sample}_Pf3D7_01_v3"
  r1="$out/${prefix}_R1.fastq.gz"
  r2="$out/${prefix}_R2.fastq.gz"
  samtools quickcheck -v "${bams[$sample]}"
  samtools view -b -P -L "$out/Pf3D7_chr1.bed" "${bams[$sample]}" |
    samtools collate -u -O -T "$out/$sample" - |
    samtools fastq -F 0x900 -1 "$r1" -2 "$r2" \
      -0 /dev/null -s /dev/null -n -
  gzip -t "$r1" "$r2"
  n1=$(gzip -cd "$r1" | awk 'END {if (NR % 4) exit 1; print NR/4}')
  n2=$(gzip -cd "$r2" | awk 'END {if (NR % 4) exit 1; print NR/4}')
  [[ "$n1" -gt 0 && "$n1" -eq "$n2" ]]
done
```

`samtools view -P -L` selects alignments overlapping the BED and fetches mates.
`collate` groups read names before conversion; `fastq -F 0x900` excludes
secondary and supplementary alignments, discards unpaired/other reads via
`-s /dev/null -0 /dev/null`, and preserves read names with `-n`.
The original whole-genome BAMs and FASTA are not included here.

## Validation of the copied dataset

All six gzip streams passed integrity checks. FASTQ records were checked for
four-line structure, header/separator markers, sequence/quality length equality,
ACGTN sequence characters, printable Phred+33 quality characters, nonzero equal
mate counts, and corresponding read names in the same order. SHA-256 checksums
matched each copied FASTQ and BED to its source. The BED contig and interval
bounds were checked against the complete chromosome reference; the reference
matched the source byte-for-byte.
