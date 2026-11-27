# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

Test data to be used for automated testing with the nf-core pipelines

## Test data for nf-core/mag

This branch contains test data for the [nf-core/mag](https://github.com/nf-core/mag) pipeline.

## CI test data

The `test_minigut` test datasets are synthetic reads generated from a single genome of a mammalian-gut-residing Bacteroides fragilis YCH46 bacterium (NCBI accession: NC_006347).

### Synthetic dataset (v6)

`test_data/v6/` is a fully synthetic metagenome of five complete bacterial genomes from NCBI RefSeq, chosen for small size and distinct GC content, plus a 500 kb region of human chromosome 21 (`host_reference/genome.hg38.chr21_10000bp_region.fa`) at low coverage as host contamination.

| Genome       | Accession       | Organism                                        | Size (Mb) | GC (%) |
| ------------ | --------------- | ----------------------------------------------- | --------- | ------ |
| `aquifex`    | GCF_000008625.1 | Aquifex aeolicus VF5                            | 1.59      | 43.5   |
| `chlamydia`  | GCF_000008725.1 | Chlamydia trachomatis D/UW-3/CX                 | 1.04      | 41.5   |
| `tropheryma` | GCF_000007485.1 | Tropheryma whipplei str. Twist                  | 0.93      | 46.5   |
| `treponema`  | GCF_000008605.1 | Treponema pallidum subsp. pallidum str. Nichols | 1.14      | 53.0   |
| `rickettsia` | GCF_000195735.1 | Rickettsia prowazekii str. Madrid E             | 1.11      | 29.0   |

There are two samples with different abundance profiles, so binners that use differential coverage have a signal to work with.
Illumina reads (2x126 bp) are simulated with [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq) 2.0.1 (HiSeq error model) and Nanopore reads with [readfaker](https://github.com/dialvarezs/readfaker) 0.4.0, using read lengths and qualities of R10.4.1 reads basecalled with dorado 2.1.2 (`sup@v5.2.0`).
Coverage is kept low (fragmented assemblies) to keep the tests fast; the goal is to exercise the pipeline flow, not to get high-quality bins.
The scripts and inputs that generate it are in [dialvarezs/mag-v6-dataset](https://github.com/dialvarezs/mag-v6-dataset).

| File                               | Description                                                 |
| ---------------------------------- | ----------------------------------------------------------- |
| `sample{1,2}_R{1,2}.fastq.gz`      | Illumina paired-end reads                                   |
| `sample{1,2}_se.fastq.gz`          | Illumina single-end reads, at twice the paired-end coverage |
| `sample1_run{1,2}_R{1,2}.fastq.gz` | The two halves of sample1, as two sequencing runs           |
| `sample{1,2}_ont.fastq.gz`         | Oxford Nanopore reads                                       |
| `record_to_genome.tsv`             | Ground truth: reference record to genome                    |
| `md5sums.txt`                      | Checksums of the reads                                      |

Read names carry their origin: Illumina reads are named `<genome>_<record>_<n>_<m>/<mate>` by InSilicoSeq, and Nanopore reads have an `XO:Z:<genome>_<record>:<start>-<end>:<strand>` comment written by readfaker.

The samplesheets are in `samplesheets/v6/`.
`assemblies/v6/` has a MEGAHIT assembly of sample1 and a metaSPAdes assembly of sample2 for assembly input tests, made by nf-core/mag from these reads and keeping only contigs of at least 1000 bp.
`databases/cat/synthetic_cat.tar.gz` is a CAT*pack 6.0.1 database with the proteins of the five genomes plus those of \_Nanoarchaeum equitans* (GCF_000008085.1) and phage phiX174 (GCF_000819615.1), with the NCBI taxonomy cut down to their lineages.
The archaeon and the phage are not in the reads: they give root and "cellular organisms" several children, as in a real database, and without them CAT_pack 6.0.1 crashes when it adds suggestive asterisks to a lineage.

The read error correction step of SPAdes (BayesHammer, 4.1.0 and 4.3.0) stalls for hours on some inputs of this data, so tests should run SPAdes with `--spades_options '--only-assembler'`.

## Full-size test data

The `manifext.full.tsv` links to gut metagenome data of antibiotic-treated patients originating from [Bertrand et al. _Nature Biotechnology_ (2019)](https://doi.org/10.1038/s41587-019-0191-2).

| SAMPLE    | ILLUMINA READS: ENA ID | ONT READs: ENA ID |
| --------- | ---------------------- | ----------------- |
| CAPES S11 | ERR3201918             | ERR3201942        |
| CAPES S21 | ERR3201928             | ERR3201952        |
| CAPES S7  | ERR3201914             | ERR3201938        |

## Downloading test data

Due the large number of large files in this repository for each pipeline, we highly recommend cloning only the branches you would use.

```bash
git clone <url> --single-branch --branch <pipeline/modules/branch_name>
```

To subsequently clone other branches[^1]

```bash
git remote set-branches --add origin [remote-branch]
git fetch
```

## samplesheets

-`samplesheets/samplesheet.long_read.csv`: Sample sheet for testing running mag with only long read input

- `samplesheets/bins_samplesheet.csv`: Sample sheet for testing `--bins_input`, listing the four bins in `bins/`

The samplesheet with the suffix `v4` (e.g., `samplesheet.long_read.v4.csv`) is intended for the v4.0.0 release of mag. In this version, samplesheets must include the `short_reads_platform` or `long_reads_platform` fields to specify the sequencing platform for short and long reads, respectively.

## Bins

`bins/` contains four MaxBin2 bins (`*.fa.gz`) used to test starting the pipeline from pre-computed bins (`--bins_input`, profile `test_bins_input`).
They were copied from `GenomeBinning/MaxBin2/bins/` of local nf-core/mag runs on the minigut reads (`samplesheets/samplesheet.multirun.v4.csv`), so the file names are the pipeline's `<assembler>-MaxBin2-<assembly id>.<n>.fa.gz`.
The two `MEGAHIT-MaxBin2-minigut.*` bins come from a MEGAHIT assembly named after the group (`minigut`), the two `SPAdes-MaxBin2-*` bins from per-sample SPAdes assemblies (`test_minigut` and `test_minigut_sample2`).
The exact pipeline version and parameters were not recorded, and the bins don't need to be regenerated: any small set of FASTA files will do.
Two of the bins, `MEGAHIT-MaxBin2-minigut.002` and `SPAdes-MaxBin2-test_minigut.002`, are near-identical and cluster together under dereplication.
`samplesheets/bins_samplesheet.csv` puts all bins in one group (`minigut`); the group only becomes the metadata id in the pipeline and has no bearing on the result.

## Databases

We have also generated tiny versions of databases for tools that normally require very large databases

### CAT_pack

For CAT of CAT_pack we did the following, and using CAT_pack (v6.0):

```bash
## Download the FASTA of coding sequences FASTA AA of B. fragilis
curl "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta_cds_aa&id=1992822979&extrafeat=null&conwithfeat=on&hide-cdd=on&ncbi_phid=CE8C15326D6BB8C10000000006490560" -o sequence.txt

## use my scripts to filter NCBI nodes/names/acc2taxid files to just a given taxid (basically fancy iterative greps)
bash ~/bin/taxdmp_filter.sh 817 ## for nodes/naames
bash ~/bin/accession2taxid_filter.sh 817 ## for accessions

## Repair headers of B. fragilis protein FASTA so the NCBI accession can be read by CAT properly
sed 's/lcl|//g;s/_/ /2' sequence.txt > sequence_fixedheaders.txt

## Generate database using CAT_pack prepare
CAT_pack prepare --db_fasta input_files/sequence_fixedheaders.txt --names input_files/names_reduced.dmp --nodes  input_files/nodes_reduced.dmp --acc2tax input_files/accession2taxid_reduced.dmp --db_dir test/

## Test using uncompressed contigs from metaspades assembly (note: --no_stars was required for some reason but seesms to only occur when we have a single genome in there possible..)
CAT_pack contigs -c SPAdes-test_minigut_sample2.scaffolds.fa -d ../../../cat_fakedb/test/db/ -t ../../../cat_fakedb/test2/tax/ --no_stars --force
```

Which resulted in a last line of the log as

```console
[2024-11-28 16:05:33] CAT is done! 658/5,223 contigs (12.60%) have taxonomy assigned.
```

And then the `test/` directory was tarred to create `minigut_cat.tar.gz`.

### geNomad

`databases/genomad/genomad_db_subset_v1.9.tar.gz` is a reduced copy of the [official geNomad database v1.9](https://zenodo.org/records/14886553), the version `genomad download-database` fetches for geNomad 1.11.x. It keeps the same layout and file formats as the official database and contains 1,000 of the 227,918 marker profiles: 5.9 MB compressed, 14 MB unpacked. Marker-based scores and taxonomy are only meaningful for the sequences it was built for.

The 1,000 markers are the ones matched by a small set of reference sequences: the _Bacteroides fragilis_ YCH46 genome (`NC_006347.1`) and short-read assemblies derived from it, plus the _Escherichia_ phage lambda genome (`NC_001416.1`) as a positive control for virus calling. The hits were collected from the `marker` column of `genomad annotate` runs against the full database, padding with a deterministic sample (seed 42) up to 1,000. Both complete genomes and fragmented assemblies were used because partial proteins match markers that full-length ones do not.

Everything else was either subset accordingly (metadata table, taxonomy mapping, minimal marker set) or copied verbatim (ICTV taxdump, integrase database, hallmark lists):

```bash
## marker names -> MMseqs2 keys (MMseqs2 18.8cc5c)
awk -F'\t' 'NR==FNR{m[$1];next} ($2 in m){print $1}' markers.txt genomad_db/genomad_db.lookup > keys.txt
mmseqs createsubdb keys.txt genomad_db/genomad_db   out/genomad_db
mmseqs createsubdb keys.txt genomad_db/genomad_db_h out/genomad_db_h

## createsubdb symlinks .lookup/.source/_mapping/_taxonomy back to the source database,
## so those are replaced by subset files; genomad_mini_db is kept as a relative symlink
## to genomad_db with its own filtered index, exactly as in the official database

tar czf genomad_db_subset_v1.9.tar.gz -C out genomad_db
```

### GTDB-Tk

We have uploaded a copy of the official [GTDB-Tk mock database](https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/mockup_db/) to this repository for use while testing mag and the nf-core GTDB modules - this significantly improves the run-time of tests as the GTDB server can be very slow.

The database is available as a gzipped tarball at:

- (r232, latest) `databases/gtdbtk/gtdbtk_mockup_r232_20260507.tar.gz`
- (r226) `databases/gtdbtk/gtdbtk_mockup_20250422.tar.gz`

> [!NOTE]
> The r232 version of the GTDB-Tk mock database has a small patch applied, because the official one contains a mismatch between the skani index and the taxonomy files, which causes GTDB-Tk to sometimes fail. The genome `GCF_003697165.2` (_E. coli_) is present in the skani index but absent from the taxonomy files, so the following line was added to `taxonomy/bac120_taxonomy_r232_reps.tsv` and `taxonomy/gtdb_taxonomy.tsv`:
>
> ```
> RS_GCF_003697165.2    d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli
> ```
>
> For more details, see [this issue](https://github.com/Ecogenomics/GTDBTk/issues/710).

### GUNC

GUNC database is provided by the GUNC team, as a minimal test set for CI/CD pipelines.

This was introduces in GUNC v1.1.0 and can by downloaded by running `gunc download_db --db test_data`. It is also available in [Zenodo](https://zenodo.org/records/19631420).

## Broken samplesheets

For testing input validation, the `samplesheets` directory contains the `broken/` subdirectory containing samplesheets with errors that should be caught by the pipeline.

-`samplesheets/broken/assembly_samplesheet_invalid_assembler.csv`: has a mistyped and missing (mandatory) assembler entry
-`samplesheets/broken/assembly_samplesheet_invalid_id.csv`: incorrectly has a space in the sample ID name
-`samplesheets/broken/assembly_samplesheet_missing_group.csv`: missing the mandatory group column
-`samplesheets/broken/assembly_samplesheet_missing_id.csv`: missing the mandatory id column
-`samplesheets/broken/assembly_samplesheet_nonunique_fasta.csv`: has duplicate FASTA files in two rows
-`samplesheets/broken/samplesheet_empty_group.csv`: group column has a header but empty entries (must have a group)
-`samplesheets/broken/samplesheet_empty_run.csv`: group column has a header but empty entries (if column is present, can't be empty)
-`samplesheets/broken/samplesheet_lr_without_sr2.csv`: long-reads given, but missing read 2 required for hybrid assembly
-`samplesheets/broken/samplesheet_missing_r1.csv`: missing mandatory read 1 column
-`samplesheets/broken/samplesheet_missing_sample.csv`: missing mandatory sample column
-`samplesheets/broken/samplesheet_nonunique_sample_run_combination.csv`: has invalid duplicate sample-run combinations
-`samplesheets/broken/samplesheet_spaces_in_name.csv`: incorrect sample name with spaces

## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).

[^1]: From [stackoverflow](https://stackoverflow.com/a/60846265/11502856)
