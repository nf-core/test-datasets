# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the nf-core pipelines

## Introduction

nf-core is a collection of high quality Nextflow pipelines.

## Documentation
nf-core/test-datasets comes with documentation in the `docs/` directory:

01. [Add a new  test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/ADD_NEW_DATA.md)
02. [Use an existing test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/USE_EXISTING_DATA.md)

## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).

## Datasets for nf-core/magmap

### Data added to implement genome preference

This is data created to deal with [issue #187](https://github.com/nf-core/magmap/issues/187).
Used in the test config `test_sourmash_mix_dupl_species`.

Seven archaeal species selected from the Sulfolobales order, each having two genomes available at GTDB (R10-RS226):

d__Archaea;p__Thermoproteota;c__Thermoprotei_A;o__Sulfolobales;f__Acidilobaceae;g__Aeropyrum;s__Aeropyrum pernix
d__Archaea;p__Thermoproteota;c__Thermoprotei_A;o__Sulfolobales;f__Ignicoccaceae;g__Ignicoccus;s__Ignicoccus hospitalis
d__Archaea;p__Thermoproteota;c__Thermoprotei_A;o__Sulfolobales;f__Ignisphaeraceae;g__Ignisphaera;s__Ignisphaera cupida
d__Archaea;p__Thermoproteota;c__Thermoprotei_A;o__Sulfolobales;f__Pyrodictiaceae;g__Pyrodictium;s__Pyrodictium delaneyi
d__Archaea;p__Thermoproteota;c__Thermoprotei_A;o__Sulfolobales;f__Sulfolobaceae;g__Metallosphaera;s__Metallosphaera hakonensis
d__Archaea;p__Thermoproteota;c__Thermoprotei_A;o__Sulfolobales;f__Sulfolobaceae;g__Metallosphaera;s__Metallosphaera javensis
d__Archaea;p__Thermoproteota;c__Thermoprotei_A;o__Sulfolobales;f__Sulfolobaceae;g__Saccharolobus;s__Saccharolobus caldissimus

Contigs for all 14 genomes (both for each species) were downloaded from NCBI and were used to:

* Generate 1000 read pairs per genome

* The 14 read pair files were concatenated 5, 5 and 4 into sample files:

testdata/archaeal_duplicates00_1.fastq.gz
testdata/archaeal_duplicates00_2.fastq.gz
testdata/archaeal_duplicates01_1.fastq.gz
testdata/archaeal_duplicates01_2.fastq.gz
testdata/archaeal_duplicates02_1.fastq.gz
testdata/archaeal_duplicates02_2.fastq.gz

* The three read pair files where included in a new samplesheet:

samplesheets/archaeal_duplicate_genomes_per_species.csv

* The species-representative genomes were used to generate a Sourmash index (archaeal_duplicates.index.sbt.zip)

testdata/archaeal_duplicates.index.sbt.zip

* A genomeinfo csv (`archaeal_duplicates.genomes.csv`) -- that also includes a Prokka-generated GFF file -- plus CheckM, CheckM2 and GTDB-Tk files were created for the non species-representative from the official GTDB metadata. Note: Only six of the genomes were selected for `archaeal_duplicates.genomes.csv` to make sure species preference tests can be done.

testdata/archaeal_duplicates.genomes.csv
testdata/archaeal_duplicates.checkm2.tsv
testdata/archaeal_duplicates.checkm.tsv

Note: `testdata/archaeal_duplicates.gtdbtk.tsv` (all seven species) is still used by nf-core/magmap's `test_bakta` profile via `bakta_test.genomes.csv`, which references all six `local_*` genomes from this dataset directly.
The rest of this seven-species dataset was retired in favour of the smaller three-species one below (nf-core/magmap#246), which is what `test_sourmash_genome_selection`/`test_species_preference` now use.

### Reduced three-species dataset for genome-selection/species-preference tests

Created for [nf-core/magmap#246](https://github.com/nf-core/magmap/issues/246): the seven-species dataset above made `sourmash_genome_selection.nf.test` and `species_preference.nf.test` extremely slow, since each of the 7-8 selected genomes per test needed a full real Prokka annotation (no pre-computed GFF is possible for genomes fetched live from NCBI).
Reduces to 3 of the same 7 species -- reusing the same genomes, GTDB metadata, and GFFs already in this dataset -- chosen to preserve every behaviour the original dataset exercised:

* **Metallosphaera javensis** -- the only species with no local/user-provided duplicate at all (`local_GCA_021654415.1` was already excluded from the genomeinfo csv in the original dataset), so its public representative (`GCF_022064045.1`) must always be picked. Tests the "no local option available" fallback path.
* **Aeropyrum pernix** -- local genome (`local_GCF_004323575.1`, CheckM completeness 98.42) beats its public representative (`GCF_000011125.1`, GTDB completeness 97.78) on raw completeness, so `--species_preference local` and `--species_preference completeness` both keep the local genome, but `--species_preference gtdb` picks the public one anyway. This is the only species in the original 7 where `completeness` and `gtdb` preference actually disagree, so it's essential for telling those two modes apart in the test.
* **Ignisphaera cupida** -- local genome (`local_GCA_023269755.1`, CheckM2 completeness 72.98) loses to its public representative (`GCF_030186535.1`, GTDB completeness 98.1) under both `completeness` and `gtdb` preference, cleanly separating `local` mode from the other two. Also the smallest genomes of the original seven, keeping Prokka runtime down.

Unlike the original dataset, whose three read samples were built from entirely disjoint sets of species, the three new samples deliberately share species pairwise, so that `--genomeset_mode joint` vs `sample` is exercised against a genome that legitimately appears in more than one sample:

* `arc00`: Metallosphaera javensis + Aeropyrum pernix
* `arc01`: Aeropyrum pernix + Ignisphaera cupida
* `arc02`: Metallosphaera javensis + Ignisphaera cupida

For each species, 3000 150bp read pairs were simulated with `wgsim` (fixed `-S` seeds for reproducibility) from one genome per species -- the local genome's FASTA for Aeropyrum pernix and Ignisphaera cupida, and the (otherwise unused) non-representative genome FASTA `GCA_021654415.1_MjAS7_1.0_genomic.fna.gz` already in this dataset for Metallosphaera javensis, since it has no local genome -- then concatenated per sample per the pairing above:

testdata/archaeal_trio00_1.fastq.gz
testdata/archaeal_trio00_2.fastq.gz
testdata/archaeal_trio01_1.fastq.gz
testdata/archaeal_trio01_2.fastq.gz
testdata/archaeal_trio02_1.fastq.gz
testdata/archaeal_trio02_2.fastq.gz

samplesheets/archaeal_trio_genomes_per_species.csv

The three species-representative genomes' existing Sourmash signatures were extracted from `archaeal_duplicates.index.sbt.zip` and reindexed into a new, smaller index (no need to re-sketch, since the genomes themselves didn't change):

testdata/archaeal_trio.index.sbt.zip

A trimmed genomeinfo csv, CheckM, CheckM2 and GTDB-Tk file were created covering just the two species with a local genome (Aeropyrum pernix, Ignisphaera cupida), reusing the exact same rows/GFFs as the original dataset:

testdata/archaeal_trio.genomes.csv
testdata/archaeal_trio.checkm.tsv
testdata/archaeal_trio.checkm2.tsv
testdata/archaeal_trio.gtdbtk.tsv

All matches were verified with real `sourmash sketch`/`prefetch` runs (not just assumed from the original dataset's numbers) before committing: every sample's simulated reads cleanly match both their public-representative and, where applicable, local-genome signatures, well above the pipeline's Sourmash Gather threshold.
