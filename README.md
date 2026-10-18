# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

Test data to be used for automated testing with the nf-core pipelines

## Test data for nf-core/mag

This branch contains test data for the [nf-core/mag](https://github.com/nf-core/mag) pipeline.

## CI test data

The `test_minigut` test datasets are synthetic reads generated from a single genome of a mammalian-gut-residing Bacteroides fragilis YCH46 bacterium (NCBI accession: NC_006347).

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

### Bakta

`databases/bakta/bakta_db_light_subset_v6.0.tar.gz` is a reduced copy of the [official Bakta light database v6.0](https://doi.org/10.5281/zenodo.14916843), the version `bakta_db download --type light` fetches for Bakta 1.12.x. It keeps the same layout and file formats as the official database (top-level directory `db-light/`): 17 MB compressed, 38 MB unpacked, instead of 1.3 GB compressed. Annotations are only meaningful for the sequences it was built for.

The subset keeps the database entries hit by the four assemblies in `assemblies/` and the four bins in `bins/`. Bakta (v1.12.0, `--meta --debug`) was run on each of them against the full light database, and the hits were collected from the raw search tables and logs Bakta keeps in debug mode. Each component was then filtered to those hits and re-indexed with the same tools Bakta uses:

- `pscc.dmnd` and `expert-protein-sequences.dmnd`: hit sequences only (`diamond getseq` + `diamond makedb`, DIAMOND 2.2.0); `bakta.db` keeps the matching `pscc` rows (also used by sORF hits)
- `pfam` and `antifam`: hit HMMs only (`hmmpress`); `ncRNA-genes` and `ncRNA-regions`: hit CMs only (`cmconvert` + `cmfetch` + `cmpress`)
- `oric.fna`: the first 10 sequences, as no assembly has an oriC hit
- `amrfinderplus-db/`: rebuilt from the AMRFinderPlus database `2026-08-07.1` (format 4.2.0, AMRFinderPlus 4.2.7) keeping only the reference proteins/CDS and HMMs of the four genes found (`tet(Q)`, `sul2`, `emrE`, `aac(3)-II`), then `amrfinder_index`
- everything else (`sorf.dmnd`, `rRNA`, `orit.fna`, `rfam-go.tsv`, `version.json`) copied verbatim

On the same inputs, the annotations match the full light database, apart from E-values and a few extra low-scoring ncRNA hits, as Infernal E-values scale with the number of models searched. The script used is `databases/bakta/build_bakta_db_subset.py`:

```bash
bakta --db db-light --meta --debug --tmp-dir runs/<name>.tmp --output runs/<name> --prefix <name> <fasta>  # for each input, after `mkdir -p runs/<name>.tmp`
python build_bakta_db_subset.py db-light db-light-subset runs/<name> [runs/<name> ...]
```

> [!NOTE]
> The AMRFinderPlus database must match the AMRFinderPlus version shipped in the Bakta container. If the Bakta module is updated to an AMRFinderPlus release with a different database format, rebuild the subset.

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
