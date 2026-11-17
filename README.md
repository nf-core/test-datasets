# nf-core/test-datasets: `abotyper`

Test data for [nf-core/abotyper](https://github.com/nf-core/abotyper), an ABO blood group genotyping pipeline for
Oxford Nanopore amplicon reads. The pipeline's `test` and `test_full` profiles and its nf-test tests read these files.

## At a glance

| | |
| --- | --- |
| **Data** | Simulated Oxford Nanopore reads of the ABO exon 2-7 amplicon, one sample per common ABO genotype |
| **Human data** | None. Every read is generated from public ISBT allele definitions on the NG_006669.2 reference |
| **Truth** | Known by construction: `simulation/truth.tsv` lists each sample's two alleles and the reads drawn from each |
| **Reproducible** | `simulation/simulate_abo_reads.py` with a fixed seed regenerates the FASTQs byte for byte |
| **Pipeline outputs** | `per_sample_processing/` comes from nf-core/abotyper run on these reads, under the real output names |
| **Licence** | [MIT](https://opensource.org/license/mit), the licence of nf-core/abotyper |

## Why simulated reads

These reads test that the pipeline runs and types correctly. They are not a measure of its accuracy, which was
assessed on laboratory data that is not part of this branch.

The first version of this branch held small subsamples of reads from cell lines with known ABO genotypes, which a
commercial partner supplied for developing the assay. The partner asked for synthetic reads to be used for the
pipeline instead, so those files were replaced by the simulated reads described here.

Simulation also has advantages of its own:

- every genotype is exact;
- any allele in the ISBT table can be added;
- the data carry no third-party terms.

## Contents

```
abotyper/
├── ABO_<genotype>_barcode<NN>.fastq.gz   simulated reads, one sample per common ABO genotype (10 files)
├── samplesheet_test_abotyper.csv          samplesheet of the 10 samples (-profile test and test_full)
├── simulation/                            how the reads were made
│   ├── simulate_abo_reads.py
│   ├── alleles.fasta                      the allele sequences the reads were drawn from
│   └── truth.tsv                          alleles, genotype and reads per allele of each sample
├── per_sample_processing/                 pipeline outputs of two samples, input to the step-level tests
│   ├── ABO_AB_barcode02/
│   └── ABO_BO1_barcode06/
└── fixtures/                              small files for individual tests (see below)
```

## Reference

Reads are aligned to the ABO RefSeqGene [NG_006669.2](https://www.ncbi.nlm.nih.gov/nuccore/NG_006669.2), which the
ISBT recommends for ABO allele annotation. The pipeline ships this sequence as `assets/refs/NG_006669v2.fasta`.
Clair3 is restricted to the exon 2-7 region `NG_006669.2:17950-24700` (`assets/refs/abo_clair3_region.bed`).

## Simulated reads

| Sample | Allele 1 | Allele 2 | Genotype | Reads |
| --- | --- | --- | --- | --- |
| ABO_AA_barcode01 | ABO\*A1.01 | ABO\*A1.01 | AA | 100 |
| ABO_AB_barcode02 | ABO\*A1.01 | ABO\*B.01 | AB | 200 |
| ABO_AO1_barcode03 | ABO\*A1.01 | ABO\*O.01.01 | AO | 100 |
| ABO_AO2_barcode04 | ABO\*A1.01 | ABO\*O.02.01 | AO | 100 |
| ABO_BB_barcode05 | ABO\*B.01 | ABO\*B.01 | BB | 100 |
| ABO_BO1_barcode06 | ABO\*B.01 | ABO\*O.01.01 | BO | 100 |
| ABO_BO2_barcode07 | ABO\*B.01 | ABO\*O.02.01 | BO | 100 |
| ABO_O1O1_barcode08 | ABO\*O.01.01 | ABO\*O.01.01 | OO | 100 |
| ABO_O1O2_barcode09 | ABO\*O.01.01 | ABO\*O.02.01 | OO | 100 |
| ABO_O2O2_barcode10 | ABO\*O.02.01 | ABO\*O.02.01 | OO | 100 |

`simulation/truth.tsv` gives the number of reads drawn from each allele.

### How the reads were made

`simulation/simulate_abo_reads.py` (Python 3, standard library and
[abotyper-core](https://pypi.org/project/abotyper-core/0.4.1/) 0.4.1) does the following:

1. **Builds each allele.** It cuts the exon 2-7 amplicon, `NG_006669.2:17922-24675` (6,754 bp), out of NG_006669.2.
   This is the span covered by full-length reads of the laboratory assay. It then applies the coding changes of the
   allele from the ISBT ABO allele table bundled with abotyper-core 0.4.1 (ISBT release 21).
   - abotyper-core is used only to place the ISBT cDNA changes on NG_006669.2, with the exon offsets the pipeline
     uses.
   - Introns stay as in NG_006669.2.
   - The resulting sequences are in `simulation/alleles.fasta`.
2. **Draws the reads.** Each read comes from one of the sample's two alleles with probability 0.5. Half of the reads
   are reverse-complemented, and up to 25 bases are trimmed from each end.
3. **Adds nanopore-like errors**, chosen to resemble R10.4.1 super-accuracy amplicon reads (about 99% read accuracy):
   - 0.4% substitutions, 0.2% insertions and 0.3% deletions per base;
   - the length of homopolymer runs of 4 or more bases changes by -1 in 6% of reads and by +1 in 2%;
   - quality scores are around Q25, and Q3-Q10 at erroneous bases.
4. **Makes the run reproducible.** The random seed (20261117) is fixed and gzip writes no timestamp, so the files are
   byte-identical between runs. Read IDs are UUIDs drawn from the same seed.

To regenerate the reads:

```bash
pip install abotyper-core==0.4.1
python simulation/simulate_abo_reads.py --ref NG_006669v2.fasta --amplicon 17922-24675 --outdir abotyper/
```

Every sample types as its intended genotype with the pipeline (see `per_sample_processing/`).

The simulated errors approximate nanopore behaviour rather than reproduce it exactly.

Simulated alleles carry no variation beyond their ISBT changes. As a result, three samples have no heterozygous SNPs
for Clair3 to phase, and Clair3 writes no gVCF or phased VCF for them:

- `ABO_AA_barcode01`: A1.01 is identical to NG_006669.2 over the amplicon;
- `ABO_AO1_barcode03` and `ABO_O1O1_barcode08`: they differ from NG_006669.2 only by c.261del.

The pipeline then types these samples from the BAM pileup, correctly, and marks them `REVIEW`. The step-level tests
therefore use `ABO_AB_barcode02` and `ABO_BO1_barcode06`, for which Clair3 writes every file.

## `per_sample_processing/`

These are the published outputs of nf-core/abotyper for `ABO_AB_barcode02` and `ABO_BO1_barcode06`, under their
real names and in the layout the pipeline writes to `<outdir>/per_sample_processing/<sample>/` (see the pipeline's
`docs/output.md`).

They were produced by running the pipeline on all ten simulated samples:

```bash
nextflow run nf-core/abotyper -r bbf404a -profile test,docker --input samplesheet.csv --outdir results
```

- **Pipeline:** nf-core/abotyper commit [`bbf404a`](https://github.com/fmobegi/abotyper/commit/bbf404a) (branch `fix/review-round3`, [nf-core/abotyper#81](https://github.com/nf-core/abotyper/pull/81)).
- **Nextflow:** 26.10.6, with the Docker containers of the pipeline.

| Folder | Step (module) | Tool and version |
| --- | --- | --- |
| `alignment/` | `MINIMAP2_ALIGN` (`-x map-ont`), `SAMTOOLS_COVERAGE`, `SAMTOOLS_FLAGSTAT`, `SAMTOOLS_STATS` | minimap2 2.30-r1287, samtools 1.25 |
| `clair3/` | `CLAIR3` (`--gvcf --enable_phasing --enable_long_indel`, model `r1041_e82_400bps_sup_v500`), `BCFTOOLS_NORM`, `BCFTOOLS_VIEW`, `BCFTOOLS_QUERY` | Clair3 2.0.5, bcftools 1.25 |
| `clair2metrics/` | `ABOTYPERCORE_CLAIR2METRICS`, `ABOTYPERCORE_PILEUP2METRICS`, `ABOTYPERCORE_CLAIR2HAPLOTYPES` | abotyper-core 0.4.1 (135 panel positions) |
| `<sample>.ABOPhenotype.txt`, `<sample>.log.txt` | `ABOTYPERCORE_PREDICTPHENOTYPE` | abotyper-core 0.4.1 |
| `mpileup/<sample>.mpileup.gz` | `SAMTOOLS_MPILEUP` (`-a -B -d 0 -Q 7 -q 1`, panel BED) | samtools 1.25 |

**`mpileup/` is not a pipeline output.** The pipeline does not publish the BAM pileup. This file was copied from the
task's work directory, unchanged, because the `ABOTYPERCORE_PILEUP2METRICS` tests need it as input.

The Clair3 VCF headers record the work-directory paths of that run.

In `<sample>.Haplotypes.tsv`, rows named `synthetic_<pos>_A0`/`_A1` are the pipeline's names for the two haplotypes
of a phase set. They are not read IDs.

## `fixtures/`

| File | Used by | Contents |
| --- | --- | --- |
| `samples_two.tsv` | `ABOTYPERCORE_AGGREGATE` tests | Sample table for `abotyper-aggregate --samples`, naming the per-sample files of `ABO_AB_barcode02` and `ABO_BO1_barcode06` |
| `samples_named.tsv` | `ABOTYPERCORE_AGGREGATE` tests | The same, with a sample name (`DonorA`) and barcode (`11`) for `ABO_AB_barcode02`, to test that they are carried into the results |
| `controls_negative.tsv` | `ABOTYPERCORE_AGGREGATE` tests | The same, with `ABO_AB_barcode02` declared a negative control, to test that a negative control with reads is flagged |
| `samplesheet_controls_multirow.csv` | pipeline test `tests/default.nf.test` | Samplesheet for the second pipeline test (see below) |
| `empty.fastq.gz` | `samplesheet_controls_multirow.csv` | A FASTQ with no reads, used as a clean negative control |

The sample tables name each file as `sample_files/<file>`. This is the folder the `ABOTYPERCORE_AGGREGATE` module
stages the per-sample files into in the pipeline and in its tests. It is not on this branch, and the files come from
`per_sample_processing/`.

Deliberate features of `samplesheet_controls_multirow.csv`:

- `ABO_BB` has barcode `7` but points at `ABO_BB_barcode05.fastq.gz`. This tests that the `barcode` column, not the
  file name, gives the reported barcode.
- `ABO_AA` has three rows pointing at the same file. This tests that rows of one sample are merged, including files
  that share a name.
- `NTC` is a negative control with barcode `95` and the empty FASTQ. This tests a clean negative control with no
  reads.
- `ABO_AO1_barcode03` has no barcode. This tests that the barcode is then taken from the sample name (`3`).
