# Metagenome binning test data

A small synthetic community containing two genomes and three samples with different coverage profiles, intended for testing metagenomic binners such as COMEBin, MetaBAT2, SemiBin, and MaxBin. It is large enough for at least one bin to exceed the minimum size that these tools apply to their final output.

## Files

| File                | Size   | Description                                     |
| ------------------- | ------ | ----------------------------------------------- |
| `contigs.fasta.gz`  | 684 kB | 449 contigs of 5 kb, 2.25 Mb total              |
| `s1.sorted.bam`     | 1.4 MB | Sample 1 reads mapped to the contigs and sorted |
| `s2.sorted.bam`     | 950 kB | Sample 2 reads mapped to the contigs and sorted |
| `s3.sorted.bam`     | 1.1 MB | Sample 3 reads mapped to the contigs and sorted |
| `s*.sorted.bam.bai` | 36 kB  | BAM indexes, required by MetaBAT2 and CONCOCT   |

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

## `multisample/`

The files above share one assembly across the three samples, which is what single-sample binners expect. Co-assembly binners that work per sample, such as `SemiBin2 multi_easy_bin` and `Vamb`, instead need one assembly per sample, combined into a single FASTA whose headers carry the sample as a prefix. `multisample/` holds that form of the same community.

Everything here is derived from the files above, so there is no second simulation to keep in step: the reads come straight back out of the published BAMs, and the contigs are the published contigs stitched and subset.

It has five samples rather than three, because SemiBin2 only accepts abundance tables in place of BAMs from five samples up.

| File                | Size       | Description                                                           |
| ------------------- | ---------- | --------------------------------------------------------------------- |
| `contigs.fasta.gz`  | 1.8 MB     | The five per-sample contig sets combined, headers `S<n>C<label>_b<m>` |
| `s*.sorted.bam`     | 1.2-1.9 MB | Each sample's reads mapped to the combined contigs, sorted            |
| `s*.sorted.bam.bai` | 20 kB      | BAM indexes                                                           |
| `s*.aemb.tsv`       | 8 kB       | Per-sample abundance over the split contigs, the form SemiBin reads   |
| `abundances.tsv`    | 10 kB      | One wide abundance table over the contigs, the form Vamb reads        |

Contigs keep the `hinf` and `port` labels, prefixed with the sample they belong to. So the expected bins are still readable from the headers, and each tool can recover the per-sample subsets.

The separator is `C`, which is Vamb's default and what SemiBin takes as `-s C`. The labels are lower case, so `C` appears exactly once in every header. Using `:` instead, SemiBin's own default, works for both tools but leaves a colon in the names of the bins Vamb writes.

| Sample | `hinf` contigs | `port` contigs | Bases   |
| ------ | -------------- | -------------- | ------- |
| `S1`   | 29             | 5              | 1.33 Mb |
| `S2`   | 31             | 3              | 1.33 Mb |
| `S3`   | 28             | 6              | 1.33 Mb |
| `S4`   | 26             | 8              | 1.33 Mb |
| `S5`   | 27             | 7              | 1.36 Mb |

The five sets overlap heavily in sequence, since they come from the same contigs. That is also what real multi-sample data looks like, and it means reads map to several samples' copies at once. Both are properties the binners are expected to handle.

### How it was generated

```bash
# 1. pull the reads back out of the published BAMs, split by genome. Both mates
#    of a pair always come from the same genome, so filtering on the reference
#    name keeps the pairs intact.
for s in s1 s2 s3; do
    for g in hinf port; do
        samtools view -h "${s}.sorted.bam" \
            | awk -v g="^${g}_" '/^@/ || $3 ~ g' \
            | samtools collate -u -O - \
            | samtools fastq -1 "reads_${s}_${g}_1.fq" -2 "reads_${s}_${g}_2.fq" \
                -0 /dev/null -s /dev/null -n
    done
done

# 2. stitch each genome's 5 kb pieces into blocks of 8, then give each sample a
#    seeded 60% of the blocks, so the per-sample sets overlap without being
#    identical. Stitching keeps the contig count near what a real assembly
#    would give, which is what the binners' runtime scales with.
#    Sample i uses random.Random(100 + i), and contigs are renamed S<i>C<label>_b<m>.

# 3. map each sample's reads to the combined contigs
bwa index contigs.fasta
for S in S1 S2 S3 S4 S5; do
    bwa mem -t 8 contigs.fasta "map_${S}_1.fq" "map_${S}_2.fq" | samtools sort -o "${S}.sorted.bam" -
    samtools index "${S}.sorted.bam"
done

# 4. the same reads as abundance tables. SemiBin reads one file per sample over
#    the split contigs, Vamb one wide table over the contigs themselves.
SemiBin2 split_contigs -i contigs.fasta.gz -o split_out
for S in S1 S2 S3 S4 S5; do
    strobealign --aemb -t 8 split_contigs.fna "map_${S}_1.fq" "map_${S}_2.fq" > "${S}.aemb.tsv"
    strobealign --aemb -t 8 contigs.fasta "map_${S}_1.fq" "map_${S}_2.fq" > "wide_${S}.tsv"
done
# joined into abundances.tsv, with the header "contigname" followed by S1 to S5
```

Each sample takes one genome's reads from one published sample, which gives five different abundance ratios out of the three that exist:

| Sample | `hinf` reads from | `port` reads from | Read pairs |
| ------ | ----------------- | ----------------- | ---------- |
| S1     | s1                | s1                | 15287      |
| S2     | s2                | s2                | 10310      |
| S3     | s3                | s3                | 11225      |
| S4     | s3                | s2                | 14846      |
| S5     | s1                | s3                | 16139      |

Software: `samtools` 1.24, `bwa` 0.7.19-r1273, `strobealign` 0.18.0, `SemiBin2` 2.4.1.

### Tested with

| Tool     | Version | Input             | Result                                     |
| -------- | ------- | ----------------- | ------------------------------------------ |
| SemiBin2 | 2.4.1   | `--input-bam`     | 2 / 2 / 3 / 2 / 2 bins, one run per sample |
| SemiBin2 | 2.4.1   | `--abundance`     | 2 bins per sample                          |
| Vamb     | 5.0.4   | `--bamdir`        | 53 clusters, 163 split bins                |
| Vamb     | 5.0.4   | `--abundance_tsv` | 53 clusters, 163 split bins                |

SemiBin2 bins each sample separately, so its counts are per sample in `S1` to `S5` order, and each run takes about 70 seconds. It was run with `-s C --ml-threshold 0 --minfasta-kbs 0 --min-len 0 --random-seed 0`.

Vamb clusters the whole set at once and then splits the bins by sample. It was run with `-o C -e 5 -q 2 3 --minfasta 1 --seed 1`, the low epoch count nf-core/modules uses for testing, which leaves most contigs in their own bin. Raising it gives fewer, larger bins.
