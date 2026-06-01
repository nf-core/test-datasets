# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

Test data to be used for automated testing with the nf-core pipelines

## Test data for nf-core/mag

This branch contains test data for the [nf-core/mag](https://github.com/nf-core/mag) pipeline.

## CI test data

The `test_minigut` test datasets are synthetic reads generated from a single genome of a mammalian-gut-residing Bacteroides fragilis\* YCH46 bacterium (NCBI accession: NC_006347).

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

The samplesheet with the suffix `v4` (e.g., `samplesheet.long_read.v4.csv`) is intended for the v4.0.0 release of mag. In this version, samplesheets must include the `short_reads_platform` or `long_reads_platform` fields to specify the sequencing platform for short and long reads, respectively.

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

### GTDB-Tk

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

-`samplesheets/broken/assembly_samplesheet_invalid_assembler.csv`: has a mistyped and missing (mandatory) assembler entry -`samplesheets/broken/assembly_samplesheet_invalid_id.csv`: incorrectly has a space in the sample ID name -`samplesheets/broken/assembly_samplesheet_missing_group.csv`: missing the mandatory group column -`samplesheets/broken/assembly_samplesheet_missing_id.csv`: missing the mandatory id column -`samplesheets/broken/assembly_samplesheet_nonunique_fasta.csv`: has duplicate FASTA files in two rows -`samplesheets/broken/samplesheet_empty_group.csv`: group column has a header but empty entries (must have a group) -`samplesheets/broken/samplesheet_empty_run.csv`: group column has a header but empty entries (if column is present, can't be empty) -`samplesheets/broken/samplesheet_lr_without_sr2.csv`: long-reads given, but missing read 2 required for hybrid assembly -`samplesheets/broken/samplesheet_missing_r1.csv`: missing mandatory read 1 column -`samplesheets/broken/samplesheet_missing_sample.csv`: missing mandatory sample column -`samplesheets/broken/samplesheet_nonunique_sample_run_combination.csv`: has invalid duplicate sample-run combinations -`samplesheets/broken/samplesheet_spaces_in_name.csv`: incorrect sample name with spaces

## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).

[^1]: From [stackoverflow](https://stackoverflow.com/a/60846265/11502856)
