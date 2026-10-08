# HG002 ONT SV genotyping test data

Small BAM and a VCF of known structural variants (SVs), for testing SV genotyping tools that genotype a given set of SVs in a long-read alignment (e.g. Sniffles `--genotype-vcf`).

- **Source**: GIAB HG002 ONT reads, R10.4.1 flowcell (PAW70337), HAC basecalling (`dna_r10.4.1_e8.2_400bps_hac@v5.0.0`), aligned to GRCh38 (chr-prefixed contigs)
- **Regions**: four windows of 1-5 kb on chr2, chr7, chr9 and chr18, containing heterozygous deletions, an insertion cluster, homozygous deletions and a heterozygous inversion
- **Downsampled**: whole reads kept by read-name hash (`samtools view -M -s 42.<fraction>`) to ~12x in the windows
- **Cleaned**: base qualities set to `*`; `MM`, `ML`, `mv`, `ts`, `ns`, `HP` and `PS` tags, `@PG` lines and `UR` paths removed

## Files

| File                            | Size     | Description                                                |
| ------------------------------- | -------- | ---------------------------------------------------------- |
| `HG002_ont_sv_genotype.bam`     | ~0.77 MB | Coordinate-sorted reads of the four windows                |
| `HG002_ont_sv_genotype.bam.bai` | ~162 KB  | BAM index                                                  |
| `HG002_ont_sv_sites.vcf.gz`     | ~3 KB    | 15 SV records of HG002 in the windows (INS, DEL, INV, BND) |
| `HG002_ont_sv_sites.vcf.gz.tbi` | <1 KB    | VCF index                                                  |

## Sites VCF

SV records called on HG002 with long-read SV callers and merged. Each record carries `SVTYPE`, `END`, `SVLEN` (and `CHR2` where present) and the callers' `GT` for HG002 (some `./.`). `QUAL` is `.`, `<TRA>` records and records longer than 100 kb are not included.

## Expected genotypes

Genotyping the VCF on the BAM with Sniffles 2.8.1 (`sniffles --input <bam> --genotype-vcf <vcf> --vcf <out>`) gives 1 `0/0`, 7 `0/1` and 7 `1/1`.
