# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

Test data to be used for automated testing with the nf-core pipelines

> ⚠️ **Do not merge your test data to `master`! Each pipeline has a dedicated branch (and a special one for modules)**

## Introduction

nf-core is a collection of high quality Nextflow pipelines. This repository contains various files for CI and unit testing of nf-core pipelines and infrastructure.

The principle for nf-core test data is as small as possible, as large as necessary. Please see the [guidelines](https://nf-co.re/docs/contributing/test_data_guidelines) for more detailed information. Always ask for guidance on the [nf-core slack](https://nf-co.re/join) before adding new test data.

## Documentation

nf-core/test-datasets comes with documentation in the `docs/` directory:

1.  [Add a new test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/ADD_NEW_DATA.md)
2.  [Use an existing test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/USE_EXISTING_DATA.md)

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

## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).

[^1]: From [stackoverflow](https://stackoverflow.com/a/60846265/11502856)

## createtaxdb CI test specific information

### FASTA files

FASTA reference files used are as follows (where possible using the same as in nf-core modules test-datasets, which are all except the human mitochondrial genome):

| Name                                     | TaxID  | Accession     | GenBank/Refseq Assembly | GenBank ID         |
| ---------------------------------------- | ------ | ------------- | ----------------------- | ------------------ |
| Human Mitochondrial genome               | 9606   | NC_012920.1   | GCF_000001405           | J01415             |
| SARS-CoV-2 genome                        | 694009 | MT192765.1    | GCA_011545545.1         | MT192765           |
| Bacteroides fragilis genome              | 817    | NZ_CP069563.1 | GCF_016889925.1         | CP069563, CP069564 |
| Candidatus portiera aleyrodidarum genome | 91844  | NC_018507.1   | GCF_000292685.1         | CP003708           |
| Streptococcus agalactiae genome          | 1311   | NZ_CP019811   | GCF_002881355.1         | CP019811           |
| Haemophilus influenzae genome            | 727    | NZ_LS483480.1 | GCF_900478275.1         | LS483480           |

> [!NOTE]
> Some reference documentation for what each accession is from:
>
> - General: https://support.nlm.nih.gov/kbArticle/?pn=KA-03434
> - GenBank: https://support.nlm.nih.gov/kbArticle/?pn=KA-03436 - [two-letter alphabetical prefix][six digits][.][version number]
> - RefSeq: https://support.nlm.nih.gov/kbArticle/?pn=KA-03437 - [two-letter alphabetical prefix][ _ ][series of digits or alphanumerical characters][.][version number]
> - Assembly: https://support.nlm.nih.gov/kbArticle/?pn=KA-03451 - [ GC{A,F} ][ _ ][nine digits][.][version number]
> - Protein: https://support.nlm.nih.gov/kbArticle/?pn=KA-03389 - [three-letter alphabetical prefix][five digits][.][version number]

The nucleotide FASTA and protein FASTA files were downloaded using `ncbi-datasets-cli` with

```bash
## Human MT genome must be extracted from the whole assembly
## Don't include coding sequences as resulting file is too large
datasets download genome accession GCF_000001405 --include genome --filename NC_012920.1.zip --chromosomes 'MT'

for i in  GCA_011545545.1 GCF_016889925.1 GCF_000292685.1 GCF_002881355.1 GCF_900478275.1; do
    datasets download genome accession --include genome,protein $i --filename $i.zip
done

for i in *.zip; do
    unzip $i
    fnaname=$(basename $(find ncbi_dataset -name '*.fna'))
    mv ncbi_dataset/data/*/*.fna .
    mv ncbi_dataset/data/*/*.faa ${fnaname%.fna}.faa
    rm -r ncbi_dataset md5sum.txt README.md
done

rm *.zip
gzip -k * ## to have both gzipped and uncompressed versions
```

### taxonomy files

These are NCBI taxdump re-constructed files, where the entries only include those of the FASTA list files above (rather than the entire tax dump).

- Downloaded NCBI [`taxdmp.tar.gz`](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/) from Feb. 2024 and used taxonkit and csvkit to filter the files to use the lineages of interest

  ```bash
  taxonkit list --ids 9606,694009,817,91844,1311,727 | taxonkit filter -E species | taxonkit lineage -t |  cut -f 3 | sed -s 's/;/\n/g' > taxids.txt
  echo 1 >> taxids.txt
  mkdir subset
  csvtk grep -Ht -f 1 -P taxids.txt /<PATH>/<TO>/nodes.dmp > subset/nodes.dmp
  cat /<PATH>/<TO>/names.dmp | csvtk fix-quotes -t | csvtk grep -Ht -f 1 -P taxids.txt | csvtk del-quotes -t > subset/names.dmp
  ## To verify
  taxonkit list --ids 1 --data-dit subset/ -nr
  ```

- Downloaded [`nucl_gb.accession2taxid.gz`](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/) from Feb. 2024, and ran

  ```bash
  accession2taxid_filter.sh 'J01415|MT192765|CP069563|CP069564|CP003708|CP019811|LS483480'
  ```

  > ![WARNING]
  > This appeared not to be correct - the accessions (first column) in this file needs
  > to refer to the assembly accession, not the individual sequences, according to Kraken2 errors.
  > So I manually replaced the original CP numbers with NC*/NZ* suffixed numbers.

- Manually made `nucl2tax.map` file based on column 2/3 of the above

- We semi-manually make a 'custom' `prot2accession` file based on the RefSeq prot accessions in the `.faa`

  ```bash
  for i in *.faa; do
      echo "accession.version	taxid" > "$i"_protaccessions.txt
      grep '>' $i | cut -d '>' -f 2 | cut -d ' ' -f 1 >> "$i"_protaccessions.txt
  done

  sed -i '1!s/$/\t694009/g' GCA_011545545.1_ASM1154554v1_genomic.faa_protaccessions.txt
  sed -i '1!s/$/\t91844/g' GCF_000292685.1_ASM29268v1_genomic.faa_protaccessions.txt
  sed -i '1!s/$/\t1311/g' GCF_002881355.1_ASM288135v1_genomic.faa_protaccessions.txt
  sed -i '1!s/$/\t817/g' GCF_016889925.1_ASM1688992v1_genomic.faa_protaccessions.txt
  sed -i '1!s/$/\t727/g' GCF_900478275.1_34211_D02_genomic.faa_protaccessions.txt

  cat *_protaccessions.txt | sort | uniq > prot.accession2taxid_custom
  rm *_protaccessions.txt
  ```

  > ![WARNING]
  > This is not an official prot.accession2taxid file, and is only used for testing purposes.
  > It uses RefSeq IDs rather than GenBank/INDSC accession numbers for mapping to taxonomy IDs
  > as the downloaded FAA files do not have the latter.

- Created a MEGAN/MALT `.db` file from the `nucl.accession2taxid_custom` file using the following `sqllite` commands:

  ```bash
  cd data/taxonomy/
  sqlite3 megan-map-custom.db
  ```

  Then within `sqlite3` create the tables (WARNING: go step by step, copy paste whole block often leads to problems reading `.` commands):

  ```sql
  -- CREATE INFO
  CREATE TABLE info (id TEXT PRIMARY KEY, info_string TEXT, size INTEGER);
  INSERT INTO info VALUES ( "general", "Created 2025-05-05 10:30:00", 8 );
  INSERT INTO info VALUES ( "Taxonomy", "Source: nucl_gb.accession2taxid (custom)", 8 );

  -- VERIFY INFO
  SELECT * FROM info;

  -- CREATE MAPPINGS
  CREATE TABLE mappings (Accession TEXT PRIMARY KEY NOT NULL, Taxonomy INTEGER);
  CREATE TABLE mappings_import (accession ANY PRIMARY KEY NOT NULL, accession_version ANY, taxid INTEGER, gi INTEGER);

  -- CONFIGURATION
  .mode tabs
  .import --skip 1 nucl_gb.accession2taxid mappings_import
  SELECT * FROM mappings_import;

  INSERT INTO mappings (Accession, Taxonomy) SELECT accession, taxid FROM mappings_import;

  -- VERIFY MAPPINGS
  SELECT * FROM mappings;

  -- CLEANUP
  DROP TABLE mappings_import;

  .quit
  ```

  To verify

  ```sql
  SELECT * FROM info;
  SELECT * FROM mappings;
  ```

## Broken Samplesheets

To help improve schema checking, we've taking then main `test.csv`, and added a few variants which have various errors.

Each file _should_ fail and give an error message from nf-schema.

- `samplesheets/broken/test_duplicate_id_and_path.csv`: has cells in both `id` and `fasta_aa` duplicated when they should be unique
- `samplesheets/broken/test_duplicate_id_only.csv`: has cells in only `id` duplicated, when all cells should be unique
- `samplesheets/broken/test_missing_both_paths.csv`: has a row where both required `fasta_dna` and `fasta_aa` paths are missing
- `samplesheets/broken/test_missing_required_column.csv`: missing the required `taxid` column
- `samplesheets/broken/test_non_existent_file.csv`: has a path to a `fasta_dna` filepath that doesn't exist

## Test Full Data

The test full input samplesheet is based on the 2025-05-08 release of [NCBI's RefSeq genome](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/) collection.

We first install `ncbi-datasets`

```bash
$ conda create -n ncbi-datasets -c conda-forge ncbi-datasets-cli -y
$ conda activate ncbi-datasets
$ datasets --version
datasets version: 18.0.5
```

> ![WARNING]
> The download commands utilising `datasets` are assume you are using [an NCBI API key](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/api-keys/) for increase queries!
> export NCBI_API_KEY=<PUT_YOUR_API_KEY_HERE>

Then we download the assembly summary files for the relevant taxa from the NCBI FTP site.

```bash
curl -o assembly_summary_archaea.txt  https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt
curl -o assembly_summary_bacteria.txt  https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
curl -o assembly_summary_viral.txt https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt
```

We condense the files to make them smaller (bacteria is 178MB!)

```bash
for i in archaea bacteria viral; do
  ## Make taxprofiler-like file
  awk -F'\t' '{if ($12 == "Complete Genome") {print $8 "\t" $6 "\t" $1"_"$16 "\t"} }' assembly_summary_$i.txt | sed 's/ /_/g' > samplesheet_"$i".txt;
  ## Make accession list file
    cut -f 3 samplesheet_$i.txt | cut -d '_' -f 1-2 >> accessionlist_$i.txt
done

rm assembly_summary_*.txt
```

<!-- TODO CAN'T USE $16 IS ITS BROKEN FOR SOME ENTRIES - N EED TO RELY PURELY ON THE GCF_ CODES -> UPDATE DOWNSTREAM STEPS -->

We can then use NCBI datasets to download the relevant files.

To estimate the size of each download

```bash
for i in archaea viral bacteria; do
  echo "$i will be"
  datasets download genome accession --inputfile accessionlist_$i.txt --assembly-level complete --include genome,protein --no-progressbar --filename assembly_$i.zip --preview
done
```

The total will be approximately the sum of the FASTA `all_genomic_fasta` and `prot_fasta` `size_mb` entries.

Then to execute (recommended in a screen/tmux session)

```bash
## Get 'dehydrated' versions of the files
for i in archaea viral bacteria ; do
  datasets download genome accession --inputfile accessionlist_$i.txt --assembly-level complete --include genome,protein --filename assembly_$i.zip --dehydrated
done

## Then 'hydrate' to generate the actual FASTA files
for i in archaea viral bacteria ; do
  unzip assembly_$i.zip -d assembly_$i/
  datasets rehydrate --directory assembly_$i/
done
```

Annoyingly the protein version of each fasta are named just `protein.faa`.

We can rename these to replicate the genomic FASTA name with the correct protein suffix.

```bash
for type in archaea viral bacteria; do
  while read fasta; do
    fnaname=$(basename $fasta)
    base=${fnaname%%_genomic.fna}
    mv $(dirname $fasta)/protein.faa $(dirname $fasta)/"$base"_protein.faa
  done < <(find assembly_$type/ -name '*.fna' -type f)
done
```

We can then prepare a 'cache' of all the full paths to the nucleotide FASTA files, which will be used in the samplesheet.

```bash
find ~+ -type f -name "*.fna" > allfna.txt
find ~+ -type f -name "*.faa" > allfaa.txt
```

Now we can reconstruct the useable sample sheet sheet based on the original samplesheet and searching through the list of FASTA files.

```bash
for type in archaea viral bacteria ; do
  if [[ -f samplesheet_"$type"_localpaths.csv ]]; then
    echo samplesheet_"$type"_localpaths.csv exists, exiting
    break
  else
    echo "samplesheet_$type.txt does not exist, generating"
  fi
  echo "id,taxid,fasta_dna,fasta_aa" > samplesheet_"$type"_localpaths.csv
  while read line; do
      genomename=$(echo $line | cut -d ' ' -f 1)
      taxid=$(echo $line | cut -d ' ' -f 2)
      accnum=$(echo $line | cut -d ' ' -f 3 | cut -d '_' -f 1-2)
      fna=$(grep -e "$accnum" allfna.txt)
      faa=$(grep -e "$accnum" allfaa.txt)
      if [[ -z "$fna" && -z "$faa" ]]; then
        echo "missing BOTH: $genomename $accnum";
        echo "$genomename $accnum" >> missing_fnafaa_"$type".txt
      elif [[ -z "$fna" && -n "$faa" ]]; then
        echo "missing FNA: $genomename $accnum"
        echo "$genomename $accnum" >> missing_fna_"$type".txt
      elif [[ -n "$fna" && -z "$faa" ]]; then
        echo "missing FAA: $genomename $accnum"
        echo "$genomename $accnum" >> missing_faa_"$type".txt
      else
        echo "${genomename}_${accnum},$taxid,$fna,$faa" >> samplesheet_"$type"_localpaths.csv
      fi
  done < <(cat samplesheet_$type.txt)
done
```

Note that this will allow that some genomes do not have protein sequences, so they will only be built with the nucleotide based profilers.

In some cases, the NCBI datasets appears to download genomes that do not have any available sequences.
Here it was only for taxa in the bacteria list.

CHECK MISSING, OR REMOVE. TO SEE WHICH: `cat samplesheet_bacteria_localpaths.csv | grep ',,'

```bash
for type in archaea viral bacteria; do
  for dir in assembly_"$type"/ncbi_dataset/data/G_*; do
    if [[ $(du --inodes $dir) -lt 3 ]]; then
      echo $dir >> missing_"$type"_inode.txt
    fi
  done
done
```

```bash
while read line; do
  taxon=$(echo $line | cut -d ',' -f 1)
  echo $taxon
  grep -e "$taxon" samplesheet_bacteria.txt
done < <(cat samplesheet_bacteria_localpaths.csv | grep ',,')
```

<!--

```bash
To prepare the samplesheet we can the do a string replacement to get the path on the system.


for type in archaea viral bacteria ; do
  if [[ -f samplesheet_"$type"_localpaths.csv ]]; then
    echo samplesheet_"$type"_localpaths.csv exists, exiting
    break
  else
    echo "samplesheet_$type.txt does not exist, generating"
  fi
  echo "id,taxid,fasta_dna,fasta_aa" > samplesheet_"$type"_localpaths.csv
  awk -v type="$type" -F'\t' '{split($3,acc,"_"); print $1 "," $2 "," "/mnt/archgen/microbiome_sciences/reference_databases/source/refseq/genomes/2025-05-08/assembly_"type"/ncbi_dataset/data/"acc[1]"_"acc[2]"/"$3"_genomic.fna" "," "/mnt/archgen/microbiome_sciences/reference_databases/source/refseq/genomes/2025-05-08/assembly_"type"/ncbi_dataset/data/"acc[1]"_"acc[2]"/"$3"_protein.faa"}' samplesheet_$type.txt >> samplesheet_"$type"_localpaths.csv
done
```

Some (particularly) virus genomes do not have protein sequences, so we need to remove these from the samplesheet.

```bash
for type in archaea viral bacteria ; do
  echo "Reading samplesheet_"$type"_localpaths.csv"
  while read line; do
    filepathfna=$(echo $line | cut -f 3 -d ',')
    if [[ ! -f $filepathfna ]]; then
      echo $filepathfna
      echo $filepathfna >> missing_fna_"$type".txt
    fi
    filepathfaa=$(echo $line | cut -f 4 -d ',')
    if [[ ! -f $filepathfaa ]]; then
      echo $filepathfaa
      echo $filepathfaa >> missing_fna_"$type".txt
    fi
  done < <(cat samplesheet_"$type"_localpaths.csv)
done

-->
