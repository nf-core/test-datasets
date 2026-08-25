# Metagenome binning test data

A small synthetic community containing two genomes and three samples with different coverage profiles, intended for testing metagenomic binners such as COMEBin, MetaBAT2, SemiBin, and MaxBin. It is large enough for at least one bin to exceed the minimum size that these tools apply to their final output.

## Files

| File               | Size   | Description                                              |
| ------------------ | ------ | -------------------------------------------------------- |
| `contigs.fasta.gz` | 684 kB | 449 contigs of 5 kb, 2.25 Mb total                       |
| `s1.sorted.bam`    | 1.4 MB | Sample 1 reads mapped to the contigs and sorted          |
| `s2.sorted.bam`    | 950 kB | Sample 2 reads mapped to the contigs and sorted          |
| `s3.sorted.bam`    | 1.1 MB | Sample 3 reads mapped to the contigs and sorted          |
| `s*.sorted.bam.bai`| 36 kB  | BAM indexes, required by MetaBAT2 and CONCOCT            |

## Composition

Contigs are named `<label>_c<n>`, so their source genomes—and therefore the expected bins—can be identified from the FASTA headers.

| Label  | Source genome                                                                | Bases   | Contigs |
| ------ | ---------------------------------------------------------------------------- | ------- | ------- |
| `hinf` | `genomics/prokaryotes/haemophilus_influenzae/genome/genome.fna.gz`           | 1.89 Mb | 378     |
| `port` | `genomics/prokaryotes/candidatus_portiera_aleyrodidarum/genome/genome.fasta` | 355 kb  | 71      |

Both genomes are included in full, so their single-copy marker sets remain complete enough for CheckM to provide meaningful completeness estimates and for binners to score candidate clusterings with marker genes.

The genomes have contrasting abundance patterns across the three samples, so differential coverage helps distinguish them in addition to their sequence composition:

| Sample | `hinf` depth | `port` depth |
| ------ | ------------ | ------------ |
| s1     | 1.5x         | 0.5x         |
| s2     | 0.5x         | 3x           |
| s3     | 1x           | 1x           |

## How it was generated

Each genome was split into fixed 5 kb contigs, and reads were simulated at the depths listed above before being mapped back to the combined set.

```bash
# 1. contigs: first record of each genome, cut into 5 kb pieces, headers relabelled
#    hinf -> 1,890,000 bp, port -> 355,000 bp

# 2. simulate reads for each genome and sample
#    pairs = depth * length / (2 * read_length)
for s in s1 s2 s3; do
    for g in hinf port; do
        wgsim -N "${PAIRS[${s}_${g}]}" -1 100 -2 100 -e 0.001 -r 0 -R 0 -X 0 -S "$seed" \
            "sim_${g}.fasta" tmp_1.fq tmp_2.fq
        cat tmp_1.fq >> "${s}_1.fq"; cat tmp_2.fq >> "${s}_2.fq"
    done
done

# 3. map back to the combined contigs
bwa index contigs.fasta
for s in s1 s2 s3; do
    bwa mem -t 8 contigs.fasta "${s}_1.fq" "${s}_2.fq" \
        | samtools sort -o "${s}.sorted.bam" -
done
```

Read pairs per genome per sample:

| Sample | `hinf` | `port` |
| ------ | ------ | ------ |
| s1     | 14175  | 895    |
| s2     | 4725   | 5370   |
| s3     | 9450   | 1790   |

The `wgsim` seeds start at 43 and increase by one for each genome–sample combination in the loop order shown above, making the simulated reads reproducible.

Software: `wgsim` 1.24, `bwa` 0.7.19-r1273, `samtools` 1.24.

## Tested with

Bins returned by each binner, in multi-sample mode (all three BAMs) and single-sample mode (`s1.sorted.bam` alone).

| Binner     | Version | 3 samples | 1 sample |
| ---------- | ------- | --------- | -------- |
| CONCOCT    | 1.1.0   | 32        | 27       |
| COMEBin    | 1.1.0   | 1         | 4        |
| MaxBin2    | 2.2.7   | 2         | 2        |
| MetaBAT2   | 2.17    | 2         | 1        |
| MetaBinner | 1.4.4   | 2         | 2        |
| SemiBin2   | 2.4.1   | 1         | 1        |
| Vamb       | 5.0.4   | 1         | 0        |

Vamb writes no FASTA at all unless `--minfasta` is given (the runs above used `--minfasta 200000`) and its VAE needs coverage across several samples to separate contigs, so a single sample leaves it with single-contig clusters.

