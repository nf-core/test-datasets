# HG002 ONT telomere test data

Small BAM/CRAM containing 17 real ONT telomeric reads from GIAB HG002, for testing telomere analysis tools (e.g. telogator2).

- **Source**: [GIAB 2025.01 ONT release](https://ont-open-data.s3.amazonaws.com/giab_2025.01/basecalling/sup/HG002/PAW70337/calls.sorted.bam), SUP basecalling, R10.4.1 flowcell (PAW70337)
- **Regions**: Last 10 kb of chr1 (`chr1:248946422-248956422`) and chr2 (`chr2:242183529-242193529`)
- **Downsampled**: `samtools view -s 42.12` to ~17 reads
- **All reads contain TTAGGG telomere repeats**

## Files

| File                           | Size    | Description                                               |
| ------------------------------ | ------- | --------------------------------------------------------- |
| `HG002_ont_tel_sub.bam`        | ~491 KB | BAM with original GRCh38 coordinates                      |
| `HG002_ont_tel_sub.bam.bai`    | ~236 KB | BAM index                                                 |
| `HG002_ont_tel_sub.cram`       | ~309 KB | CRAM with coordinates adjusted to mini reference          |
| `HG002_ont_tel_sub.cram.crai`  | 66 B    | CRAM index                                                |
| `HG002_ont_tel_sub_ref.fa`     | ~149 KB | Mini reference: last ~56 kb of chr1 + last ~94 kb of chr2 |
| `HG002_ont_tel_sub_ref.fa.fai` | 42 B    | Mini reference index                                      |

## CRAM mini-reference

The original BAM is aligned to GRCh38, but including full chromosome sequences would be too large for this repo. Instead, the CRAM uses a mini reference containing only the telomeric regions where the reads align:

- `chr1`: positions 248,900,001-248,956,422 of GRCh38 (56,422 bp)
- `chr2`: positions 242,100,001-242,193,529 of GRCh38 (93,529 bp)

Read coordinates in the CRAM are offset accordingly (e.g. GRCh38 chr1:248,939,253 becomes chr1:39,253 in the CRAM). The read sequences and qualities are identical to the BAM.

## How to regenerate

### BAM

Requires Docker and internet access:

```bash
IMAGE="community.wave.seqera.io/library/telseq_samtools_bamtools:0674b50fa5a173e1"
BAM="https://ont-open-data.s3.amazonaws.com/giab_2025.01/basecalling/sup/HG002/PAW70337/calls.sorted.bam"

docker run --rm -v "$(pwd):/out" "${IMAGE}" bash -c "
  samtools view -b '${BAM}' 'chr1:248946422-248956422' 'chr2:242183529-242193529' -o /tmp/tel_full.bam
  samtools index /tmp/tel_full.bam
  samtools view -s 42.12 -b /tmp/tel_full.bam -o /out/HG002_ont_tel_sub.bam
  samtools index /out/HG002_ont_tel_sub.bam
"
```

### CRAM

Requires the BAM above plus GRCh38 chr1 and chr2 FASTA sequences. See `build_ont_cram.py` in the [telomerelength pipeline repo](https://github.com/seqera-services/m42-telomerelength) for the conversion script, which:

1. Extracts telomeric regions from GRCh38 into a mini reference FASTA
2. Rewrites BAM coordinates (subtracting region offsets) using pysam
3. Converts to CRAM using the mini reference
