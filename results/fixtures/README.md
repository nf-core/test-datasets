# GWAS pipeline fixtures

Pipeline-level test fixtures for `nf-core/gwas`: one cohort, all 22 autosomes, three
genotype encodings, one signal-bearing quantitative trait, one signal-bearing binary
trait, and a covariate set with sex, age and principal components.

They exist because the older `results/pheno_cov/` files cannot exercise the current
routes: they are VCF-only, single-chromosome, headerless, have no binary trait, no
sex/age/PC covariates, and a phenotype that is uncorrelated noise, which leaves
variance-component estimators on the boundary of the parameter space.

## Provenance

Everything here is derived from `results/chunked_vcfs/chr{1..22}_chunked.vcf.gz`,
which are already committed on this branch (1000 Genomes phase 3, GRCh37/b37, 2504
samples, the first 4500 variant records of each chromosome). **No step downloads
anything.**

## Contents

### `genotypes/`

Twenty-three bundles — one per autosome plus the merge — each present in all three
encodings, describing exactly the same samples and the same variants:

| bundle          | chromosome | samples | variants |
| --------------- | ---------- | ------- | -------- |
| `example_chr1`  | 1          | 500     | 319      |
| `example_chr2`  | 2          | 500     | 172      |
| `example_chr3`  | 3          | 500     | 181      |
| `example_chr4`  | 4          | 500     | 181      |
| `example_chr5`  | 5          | 500     | 190      |
| `example_chr6`  | 6          | 500     | 269      |
| `example_chr7`  | 7          | 500     | 261      |
| `example_chr8`  | 8          | 500     | 259      |
| `example_chr9`  | 9          | 500     | 203      |
| `example_chr10` | 10         | 500     | 167      |
| `example_chr11` | 11         | 500     | 273      |
| `example_chr12` | 12         | 500     | 274      |
| `example_chr13` | 13         | 500     | 175      |
| `example_chr14` | 14         | 500     | 82       |
| `example_chr15` | 15         | 500     | 178      |
| `example_chr16` | 16         | 500     | 209      |
| `example_chr17` | 17         | 500     | 263      |
| `example_chr18` | 18         | 500     | 223      |
| `example_chr19` | 19         | 500     | 231      |
| `example_chr20` | 20         | 500     | 185      |
| `example_chr21` | 21         | 500     | 362      |
| `example_chr22` | 22         | 500     | 285      |
| `example_all`   | 1–22       | 500     | **4942** |

Encodings per bundle:

- PLINK 2: `.pgen` + `.pvar` + `.psam`
- PLINK 1: `.bed` + `.bim` + `.fam`
- VCF: `.vcf.gz` + `.vcf.gz.tbi` (bgzip-compressed, tabix-indexed, phased `GT`)

Both the per-chromosome split and the merged bundle are published: REGENIE and
LDAK-KVIK need more than one chromosome for leave-one-chromosome-out to be
non-degenerate, and the single-cohort routes want one merged bundle.

Sample identifiers are the bare 1000 Genomes IDs (`HG00096`, …) in every encoding:
`.psam` `IID`, `.fam` column 2, and the VCF `#CHROM` sample columns all carry the same
string, and it matches the `IID` column of everything under `pheno_cov/`. The PLINK
filesets also carry an `FID`, set equal to the `IID`; the VCFs are exported with
`--export vcf id-paste=iid` so the `FID` is not pasted onto the sample name. That
matters for the VCF association route, where the genotypes are passed as `--vcf` and
joined to a phenotype file by `IID`.

Sample sex is recorded in `.psam` and `.fam` (1 = male, 2 = female).

### `pheno_cov/`

| file                  | header                     | description                                             |
| --------------------- | -------------------------- | ------------------------------------------------------- |
| `example.pheno`       | `FID IID QT BT`            | quantitative trait and binary trait                     |
| `example.covar`       | `FID IID SEX AGE PC1..PC5` | all covariates in one file                              |
| `example.qcovar`      | `FID IID AGE PC1..PC5`     | quantitative covariates only                            |
| `example.catcovar`    | `FID IID SEX`              | categorical covariates only                             |
| `fixture_samples.tsv` | `#FID IID`                 | the selected 500 samples                                |
| `fixture_sex.tsv`     | `#FID IID SEX`             | the sex assignment                                      |
| `causal_weights.tsv`  | `ID A1 BETA_QT BETA_BT`    | the simulated causal effects (ALT is the effect allele) |
| `fixture_summary.txt` | —                          | realised properties of the simulation                   |

All files are tab-delimited. `BT` uses the PLINK convention, `1` = control,
`2` = case. GCTA and LDAK want headerless phenotype files, so the pipeline is
expected to normalise these before use. `gcta --fastGWA-mlm-binary` additionally
wants `0`/`1` rather than `1`/`2` and errors with `XtWX is not invertible` on the
PLINK coding, so that route needs the recode as well as the header strip.

## How the data were built

1. **Samples.** The 2504 source sample IDs are sorted with `LC_ALL=C sort` and
   every 5th is taken until 500 are selected. There is no RNG in the selection.
2. **Sex.** The source VCFs carry no sex, and the committed chrX chunk lies
   entirely in PAR1 where both sexes are diploid, so sex cannot be recovered from
   the data; it is assigned by alternating over the sorted selection, giving 250
   males and 250 females.
3. **Variants.** Per chromosome, `plink2` keeps biallelic ACGT SNPs with unique
   IDs and MAF ≥ 0.005, then LD-thins with `--indep-pairwise 100 10 0.2`. `INFO`,
   `QUAL` and `FILTER` are dropped, because the source `INFO` describes all 2504
   samples and is stale for a 500-sample subset. Variants without an rsID are
   named `CHROM:POS:REF:ALT`.
4. **Merge and PCA.** The 22 per-chromosome bundles are merged with
   `--pmerge-list`, and ten principal components are computed from the merged
   bundle with exact (non-approximate) `plink2 --pca`.
5. **Genetic values.** Two disjoint sets of 150 causal variants are drawn from the
   merged bundle with effect sizes drawn from N(0, 1) under a fixed seed, one set
   for the quantitative trait and one for the binary trait's liability. Each
   sample's genetic value is computed from the actual genotypes with
   `plink2 --score`.
6. **Traits.** Both the quantitative trait and the binary liability are built as

   ```
   sqrt(0.5) * standardised genetic score
   + sqrt(0.1) * standardised covariate effect (sex and age)
   + sqrt(0.4) * N(0, 1) noise
   ```

   so the target heritability is 0.5 on the total-variance scale and the supplied
   covariates genuinely matter. The binary trait dichotomises the liability at the
   empirical quantile that gives exactly the target prevalence, so the case count
   is fixed rather than sampled.

## Chromosome coverage and the LD threshold

The fixture covers all 22 autosomes rather than a handful, and that is what makes
the LD thinning real. The source chunks are the _first_ 4500 variant records of each
chromosome, which span 0.773 Mb on chr1 but only 0.128 Mb on chr2 and 0.104 Mb on
chr3 — 100–130 kb windows holding 4500 variants each, i.e. essentially one LD block
per chromosome. Within such a window independence is absent from the data rather
than merely under-filtered, so no r² threshold recovers it: the number of genuinely
independent variants tracks the number of _chromosomes_ used, not the number of
variants kept. Cross-chromosome variants are unlinked by construction, so widening
to 22 autosomes buys independence proportional to chromosome count while
simultaneously lifting the variant count.

A hard floor applies: `gcta/fastgwa`'s binary-trait path needs at least 2000
autosomal SNPs for tuning and aborts with `can't read 2000 SNPs from the autosomes
for tuning` below that.

Generation started at `100 10 0.2`, with relaxation required only if the 2000-variant
floor was missed:

| `--indep-pairwise` | chromosomes | surviving variants | verdict                                     |
| ------------------ | ----------- | ------------------ | ------------------------------------------- |
| `100 10 0.2`       | 1–22        | **4942**           | **chosen** — conventional r², clears 2000   |
| `100 10 0.9`       | 1–3         | 1706               | previous fixture; below the floor, rejected |

`100 10 0.2` is the conventional pruning threshold and was the initial
point; because it cleared the floor by a wide margin on that first run, no relaxation
was needed and none was applied. It ships because it is that starting point and it
worked, not because tighter thresholds were run and 0.2 came out best — no threshold
below 0.2 was tested against the 22-autosome input. As a known-headroom note rather
than a systematic sweep, `100 10 0.1` was independently run once against the same
input and also clears the floor, at 3801 variants; nothing tighter than 0.1 has been
tried. The 0.9/three-chromosome row is the fixture this one replaces, kept as the
measurement that motivated the change.

Verified against `gcta64` directly: `--fastGWA-mlm-binary` on the previous
1706-variant merged bundle fails with `can't read 2000 SNPs from the autosomes for
tuning`, while on this 4942-variant bundle it completes and writes results for all
4942 SNPs.

## Realised properties

From `pheno_cov/fixture_summary.txt`:

```
samples                       500
cases / controls              150 / 350   (prevalence 0.30)
males / females               250 / 250
cor(QT, genetic score)        0.6983   ->  r2 = 0.4876
cor(liability, genetic score) 0.7106
cor(BT, genetic score)        0.5231
mean age                      54.84
```

A `plink2 --glm` run of `QT` on `example_all` with `SEX`, `AGE` and `PC1..PC5` as
covariates gives a minimum p-value of 7.5e-6 over the 4942 tests, one variant below
the Bonferroni threshold (0.05/4942 = 1.0e-5) and ten below 1e-3. Per-variant
p-values are weaker than the previous three-chromosome fixture's 1.0e-12 precisely
because the thinning is now genuine: at r² = 0.2 the 150 causal variants are no
longer shadowed by near-duplicate tags, and the same simulated heritability is
spread over 4942 independent tests instead of 1706 correlated ones. The trait-level
signal is unchanged (the genetic score still explains 49% of `QT`'s variance), so
heritability and mixed-model routes see the same effect; only single-variant
significance is diluted.

## Regenerating

From the root of this branch:

```bash
nextflow run main.nf -profile test --skip_source_generation
```

`--skip_source_generation` skips the 1000 Genomes download and chunking stage,
whose output is already committed under `results/chunked_vcfs/`. Everything under
`results/fixtures/` is rewritten in place, as real copies rather than symlinks.

The output is byte-for-byte reproducible: every RNG draw goes through a fixed seed
with R's RNG kind pinned explicitly, `plink2` is run single-threaded with a fixed
seed, sample selection and merge order are deterministic, and the run-date header
`plink2` writes into exported VCFs is stripped before compression. Two independent
runs of the command above, into separate output and work directories, produce
byte-identical copies of all 192 generated files.

The simulation is controlled by the `fixture_*` parameters in `nextflow.config`
(chromosomes, sample count, MAF floor, LD thinning, causal-variant count,
heritability, covariate variance, prevalence, seed). Changing any of them changes
the fixtures, and any pipeline snapshot taken against them has to be regenerated.
