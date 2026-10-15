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

### NCBI assembly summaries and GTDB metadata for the three-species dataset

Created for [nf-core/magmap#258](https://github.com/nf-core/magmap/issues/258): the tests using `archaeal_trio.index.sbt.zip` read NCBI's full assembly summaries (about 2 GB) and GTDB's full archaeal metadata.
Downloading and parsing the summaries dominated their runtime, and the GTDB download has made test runs fail when the server was unavailable.

These are extracts of the real files, keeping the rows the tests need plus decoys, so that the lookups have to pick the right rows rather than the only ones:

* The two header lines of NCBI's `assembly_summary_refseq.txt` and `assembly_summary_genbank.txt` (downloaded 2026-09-23), with:
  * the three public genomes in the index (`GCF_000011125.1`, `GCF_022064045.1`, `GCF_030186535.1`), and their GenBank twins (`GCA_` with the same number);
  * up to four other genomes per genus (Aeropyrum, Metallosphaera, Ignisphaera), including viruses named after the genus;
  * two GenBank rows whose `ftp_path` is `na`;
  * every 250,000th row of each file.
* One constructed row: `GCA_977173205.1` in the GenBank file is a real row truncated before the `ftp_path` column, to reproduce the missing-field case from nf-core/magmap#244.
* The header of GTDB release 226 `ar53_metadata_r226.tsv.gz`, with every genome of the three species, up to three other genomes per genus, and every 1500th row.
  This includes the GTDB rows for the local genomes' own accessions and for the non-representative Metallosphaera javensis genome.

All rows were extracted from the live files with the commands below; none except the truncated one were edited.
The NCBI summaries change daily, so rerunning the first command picks different sampled rows.

```bash
# NCBI assembly summaries (downloaded 2026-09-23)
for src in refseq genbank; do
    curl -s https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_${src}.txt | awk -F'\t' '
        NR<=2 { print; next }
        $1 ~ /^GC[AF]_(000011125|022064045|030186535)\./ { print; next }
        $8 ~ /^(Aeropyrum|Metallosphaera|Ignisphaera) / && g[substr($8,1,index($8," "))]++ < 4 { print; next }
        ($20 == "" || $20 == "na" || NF < 20) && bad++ < 2 { print; next }
        NR % 250000 == 0 { print }
    ' > trio_${src}.txt
done
cp trio_refseq.txt archaeal_trio.assembly_summary_refseq.txt

# Truncate one GenBank row before the ftp_path column (column 20)
awk -F'\t' -v OFS='\t' '$1=="GCA_977173205.1"{NF=19} {print}' trio_genbank.txt > archaeal_trio.assembly_summary_genbank.txt

# GTDB release 226 archaeal metadata (gawk, for gensub)
curl -s https://data.gtdb.ecogenomic.org/releases/release226/226.0/ar53_metadata_r226.tsv.gz | zcat | gawk -F'\t' '
    NR==1 { print; next }
    $20 ~ /s__(Aeropyrum pernix|Metallosphaera javensis|Ignisphaera cupida)$/ { print; next }
    $20 ~ /g__(Aeropyrum|Metallosphaera|Ignisphaera);/ && g[gensub(/.*g__([^;]+);.*/, "\\1", 1, $20)]++ < 3 { print; next }
    NR % 1500 == 0 { print }
' > archaeal_trio.ar53_metadata.tsv
```

testdata/archaeal_trio.assembly_summary_refseq.txt
testdata/archaeal_trio.assembly_summary_genbank.txt
testdata/archaeal_trio.ar53_metadata.tsv

### HMM profiles for the HMM search

Used in the `-profile test, hmmfiles` test in `tests/default.nf.test`, added for [issue #34](https://github.com/nf-core/magmap/issues/34).

Three Pfam profiles downloaded from [InterPro](https://www.ebi.ac.uk/interpro/) that hit proteins of `GCA_002688505.1`, the genome in the default test:

* `PF00006.hmm`: ATP-synt_ab, ATP synthase alpha/beta family, nucleotide-binding domain (PF00006.32)
* `PF00009.hmm`: GTP_EFTU, elongation factor Tu GTP binding domain (PF00009.34)
* `PF00154.hmm`: RecA_N, recA bacterial DNA recombination protein (PF00154.28)

### Output of nf-core/mag with dereplication

Created for `--mag_outdir` ([nf-core/magmap#280](https://github.com/nf-core/magmap/issues/280)), which reads the output directory of an nf-core/mag run with `--dereplicate`.
It is a trimmed copy of such an output directory, with the same layout as the real one under `testdata/mag_dereplicated/`.

The input was four bins named after assembler, binner and sample (the `minigut` test samples of nf-core/mag), run through nf-core/mag 5.6.0dev with `--bins_input` and `--dereplicate`.
`--bins_input` was an unreleased branch of nf-core/mag at the time of writing: commit `f5218a59` on the branch `bins-input` of [erikrikarddaniel/nf-core-mag](https://github.com/erikrikarddaniel/nf-core-mag).
The run used, among others, `--bins_input bins.csv`, `--binning_map_mode group`, `--run_checkm2 --checkm2_db binette_checkm2_small_db.dmnd`, `--gtdb_db gtdbtk_mockup_r232_20260507.tar.gz`, `--skip_prodigal`, `--dereplicate --dereplicate_ani 90 --dereplicate_min_completeness 1 --dereplicate_max_contamination 100` (so that the tiny bins pass the filters) and `--gtdbtk_min_completeness 0.01 --gtdbtk_min_perc_aa 10`.
Only the files magmap reads were kept.
If the layout of the output changes when `--bins_input` is released, rerun with the released version and replace the files.

The four bins form three clusters: `MEGAHIT-MaxBin2-minigut.002` is a member of the cluster of `SPAdes-MaxBin2-test_minigut.002`.
So there are three representatives, and the fourth bin is in neither the dereplicated bins nor the GTDB-Tk and Prokka results, but is in the CheckM2 summary.
The GTDB-Tk summary covers a species-level call, a partial lineage (`g__;s__`) and one `Unclassified Bacteria`.
GTDB-Tk used a mock-up database, so the taxonomy is only meant for testing the parsing: `SPAdes-MaxBin2-test_minigut.002` is classified as _Escherichia coli_ although its closest placement is _Selenomonas_C_.

```text
testdata/mag_dereplicated/
  GenomeBinning/Dereplication/dereplicated_bins/*.fa.gz
  GenomeBinning/Dereplication/dereplicated_bins.tsv
  GenomeBinning/QC/checkm2_summary.tsv
  Taxonomy/GTDB-Tk/gtdbtk_summary.tsv
  Annotation/Prokka/external/<bin>/<bin>.gff
```

The GFF files are the Prokka output cut before the `##FASTA` section.
The sequences are already in the dereplicated bins, and magmap only reads the annotation lines.
