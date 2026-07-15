# Test data for `gedi/price`

A minimal cohort of four Ribo-seq samples covering protein-coding regions of chr19 + chr22 (Ensembl GRCh38). Genome reference and BAMs are masked/filtered to the protein-coding loci so that PRICE's expectation-maximisation converges and produces a non-empty `orfs.tsv` while keeping every fixture file under 4 MiB and the whole set under 11 MB.

## Why this minimal set?

PRICE's ORF inference fails (`Index 0 out of bounds for length 0` in `PriceOrfInference`) when too few candidate ORFs feed the noise model. Empirically:

- chr19 alone, six samples: PRICE estimates the model but returns zero ORFs.
- chr19 + chr22, three samples: PRICE crashes in `PriceOrfInference`.
- chr19 + chr22, four samples, protein-coding-only reference: 381 ORFs - works.

Restricting the reference to protein-coding genes (drops ~10% of ORF calls vs the full annotation), filtering BAMs to protein-coding gene loci (drops reads PRICE cannot use anyway), and keeping the four deepest samples is the minimum that keeps the test meaningful.

## Files

| File                                                | Size            | Description                                                                                                                                                                                           |
| --------------------------------------------------- | --------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Homo_sapiens.GRCh38_chr19_22.pc_exon_masked.fa.gz` | 3.3 MB          | chr19+chr22 from Ensembl GRCh38 primary assembly; everything outside protein-coding-gene exons hard-masked to `N`.                                                                                    |
| `Homo_sapiens.GRCh38.111_chr19_22.pc.gtf.gz`        | 1.6 MB          | chr19+chr22 from Ensembl 111 GTF, subset to `gene_biotype "protein_coding"`. Attribute column trimmed to `gene_id`, `transcript_id`, `gene_biotype`, `gene_name`, `transcript_biotype`.               |
| `bams/SRX1178088{5,6,7,8}.chr19_22.ds50.bam`        | 1.0-1.5 MB each | Four Ribo-seq samples from GSE182201, STAR-aligned to GRCh38, filtered to chr19+chr22, downsampled to 50% with `samtools view -bs`, then filtered to reads overlapping protein-coding-gene intervals. |
| `bams/*.bai`                                        | <90 KB each     | BAM indexes.                                                                                                                                                                                          |

Total: ~11 MB across 10 files.

## How they were derived

1. Source Ribo-seq reads SRR15480788/9/90/91 (SRA accessions for SRX11780885-8, from GSE182201) were aligned to the Ensembl GRCh38 primary assembly + Ensembl 111 GTF using STAR (full-genome alignment, post-rRNA-filtering).
2. Each genome-sorted BAM filtered to chr19+chr22 with `samtools view -bh -F 256 <bam> 19 22`.
3. Downsampled to 50% with `samtools view -bh -s 42.50`.
4. Filtered to reads overlapping protein-coding gene loci with `samtools view -bh -L pc_gene.bed`.
5. FASTA built from `Homo_sapiens.GRCh38.dna.chromosome.{19,22}.fa.gz` (Ensembl 111), then N-masked outside the union of protein-coding-gene exon intervals.
6. GTF subset to chr19+chr22 from `Homo_sapiens.GRCh38.111.chr.gtf.gz`, restricted to `gene_biotype "protein_coding"`, then stripped of non-essential attributes.

## Verified

Running `gedi -e IndexGenome` then `gedi -e Price` on this cohort yields a `run.orfs.tsv` with 381 lines (380 ORFs plus header).

Used by `modules/nf-core/gedi/price/tests/main.nf.test`.
