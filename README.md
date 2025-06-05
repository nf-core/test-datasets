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

The full test data here is designed to correspond to the test data used for the full test run of [nf-core/taxprofiler](https://nf-co.re/taxprofiler).

FASTA files for use in database construction are taken from Supp. Table 1 from [Meslier (2022)](https://doi.org/10.1038/s41597-022-01762-z).
The list of genomes were copy-pasted from the Excel document into an empty tsv file called `41597_2022_1762_MOESM1_ESM.tsv`.

NCBI Datasets package (v18.05) was then used to download the reference genomes and protein translations of each strain from the file.

```bash
awk -F'\t' '{print $4}' 41597_2022_1762_MOESM1_ESM.tsv | tail -n +2 > download_urls.txt

datasets download genome accession --inputfile download_urls.txt --include genome,protein --filename meslier_2022.zip --preview
```

To download

```bash
datasets download genome accession --inputfile download_urls.txt --include genome,protein --filename meslier_2022.zip --dehydrated
```

To unpack and retrieve genomes

```bash
unzip meslier_2022.zip
datasets rehydrate --directory .
```

We can validate the number of genome files downloaded with:

```bash
find -name '*.fna' -type f | wc -l ## 89, expecting 91
find -name '*.faa' -type f | wc -l ## 87, expecting 91
```

The 2 fna discordance is because we assume the assemblies have been suppressed by NCBI, and thus not available for download.
Thus we continue with what we have.

Annoyingly the protein version of each fasta are named just `protein.faa`.

We can rename these to replicate the genomic FASTA name with the correct protein suffix.

```bash
while read fasta; do
  fnaname=$(basename $fasta)
  base=${fnaname%%_genomic.fna}
  mv $(dirname $fasta)/protein.faa $(dirname $fasta)/"$base"_protein.faa
done < <(find ncbi_dataset/ -name '*.fna' -type f)
```

During this renaming, two files were reported as missing.

```bash
mv: cannot stat 'ncbi_dataset/data/GCA_000308215.1/protein.faa': No such file or directory
mv: cannot stat 'ncbi_dataset/data/GCA_002563335.1/protein.faa': No such file or directory
```

Again, not sure why.
But indeed on the NCBI FTP server ([GCA_000308215.1](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/308/215/GCA_000308215.1_ASM30821v1/), [GCA_002563335.1](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/563/335/GCA_002563335.1_ASM256333v1/)) there is no protein translation for the these GenBank assemblies (they appear to exist for the RefSeq versions, but we will stay with GenBank for consistency. Maybe submitted assemblies were not originally annotated and RefSeq did it for them.).

We can then prepare a samplesheet of all the FASTA files, which will be used as input to the pipeline.

```bash
echo "id,taxid,fasta_dna,fasta_aa" > samplesheet.csv
while read accession; do
  echo "Preparing: $accession"
  taxid=$(datasets summary genome accession $accession | jq '.reports[0].organism.tax_id')
  fna=$(find ~+ -name "*${accession}*.fna" -type f)
  faa=$(find ~+ -name "*${accession}*.faa" -type f)
  if [[ -z "$fna" && -z "$faa" ]]; then
    echo "missing BOTH: $accession"
    echo "$accession" >> missing_fnafaa.txt
  else
    echo "$accession,$taxid,$fna,$faa" >>  samplesheet.csv
  fi
done < download_urls.txt
```

We then download the relevant taxonomy files from NCBI, which are used by nf-core/createtaxdb.

```bash
## NCBI accession2taxid files
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz ## 464617.pts-70.bionc21

## NCBI taxdmp files
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar xvzf taxdump.tar.gz
```

We can reconstruct the 'non-standard' taxonomy files from the above with:

```bash
## Command borrowed from from https://github.com/khyox/recentrifuge/wiki/Centrifuge-nt#step-by-step-instructions
gunzip -c nucl_gb.accession2taxid.gz | awk -v OFS='\t' '{print $2, $3}' >> nucl2taxid.map
gunzip -c prot.accession2taxid.gz | awk -v OFS='\t' '{print $2, $3}' >> prot2taxid.map
```

We will also download the relevant MALT files:

```bash
wget https://software-ab.cs.uni-tuebingen.de/download/megan6/nucl_acc2tax-Jul2019.abin.zip
unzip nucl_acc2tax-Jul2019.abin.zip
```

Test command (to be replaced with config):

```bash
nextflow run nf-core/createtaxdb -r dev -profile mpcdf_viper --input samplesheet_viper.csv --outdir ./results --dbname test_full --accession2taxid nucl_gb.accession2taxid.gz --prot2taxid prot2taxid.map --nucl2taxid nucl2taxid.map --nodesdmp nodes.dmp --namesdmp names.dmp --malt_mapdb nucl_acc2tax-Jul2019.abin --malt_mapdb_format a2t --build_bracken --build_centrifuge --build_diamond --build_ganon --build_kaiju --build_kraken2 --build_krakenuniq --krakenuniq_build_options '--jellyfish-bin $(which jellyfish)' --build_malt --build_kmcp --generate_tar_archive --generate_pipeline_samplesheets taxprofiler --generate_samplesheet_dbtype tar --unzip_batch_size 100 -c custom.config
```

<!-- TODO MOVE TO RIGHT PLACE -->

<!-- OLD VERSION - DO NOT USE TOO BIG

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

TODO CAN'T USE $16 IS ITS BROKEN FOR SOME ENTRIES - N EED TO RELY PURELY ON THE GCF_ CODES -> UPDATE DOWNSTREAM STEPS

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

We can check these by checking the number of filesystem objects in each of the directories.
We expect 3: the directory itself, then a FNA and FAA file.
Anything less than this indicates that at least one of the files is missing.

```bash
for type in archaea viral bacteria; do
  for dir in assembly_"$type"/ncbi_dataset/data/G_*; do
    if [[ $(du --inodes $dir) -lt 3 ]]; then
      echo $dir >> missing_"$type"_inode.txt
    fi
  done
done
```

We can then re-run the above steps above in a condensed manner to re-download for missing files.

For missing FAA files, we can use the NCBI datasets API to check if the protein FASTA files actually do exist

```bash
while read line; do
  query=$(echo $line | cut -d ' ' -f 2)
  baseurl=$(curl --silent -X GET "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$query/links"  -H 'accept: application/json' | jq '.[] | .[0] |.resource_link')
  assembly=$(echo $baseurl | rev | cut -d '/' -f 1 | rev)
  fullurl=${baseurl//\"/}/${assembly//\"/}_protein.faa.gz

  if curl --output /dev/null --silent --head --fail "$fullurl"; then
    echo "valid: $fullurl" >> missing_faa_viral_validation.txt
  else
    echo "missing: $fullurl" >> missing_faa_viral_validation.txt
  fi

done < missing_faa_viral.txt

## These should be redownloaded
grep 'valid:' missing_faa_viral_validation.txt
```

In this case all of those genomes do not have protein sequences, so we can ignore them assuming we want a truly comparable dataset across all profiling databases.

We can also double check files from species where both FNA and FAA are truly missing

```bash
while read line; do
  query=$(echo $line | cut -d ' ' -f 2)
  ls -l assembly_bacteria/ncbi_dataset/data/$query/*
done < missing_fnafaa_bacteria.txt
```

For missing FNA/FAA files, we can use the NCBI datasets API to check if the genomic FASTA files actually do exist

```bash
while read line; do
  query=$(echo $line | cut -d ' ' -f 2)
  baseurl=$(curl --silent -X GET "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$query/links"  -H 'accept: application/json' | jq '.[] | .[0] |.resource_link')
  assembly=$(echo $baseurl | rev | cut -d '/' -f 1 | rev)
  fullurlfna=${baseurl//\"/}/${assembly//\"/}_genomic.fna.gz
  fullurlfaa=${baseurl//\"/}/${assembly//\"/}_protein.faa.gz

  if curl --output /dev/null --silent --head --fail "$fullurlfna"; then
    echo "valid: $fullurlfna" >> missing_fnafaa_bacteria_validation.txt
  else
    echo "missing: $fullurlfna" >> missing_fnafaa_bacteria_validation.txt
  fi

  if curl --output /dev/null --silent --head --fail "$fullurlfaa"; then
    echo "valid: $fullurlfaa" >> missing_fnafaa_bacteria_validation.txt
  else
    echo "missing: $fullurlfaa" >> missing_fnafaa_bacteria_validation.txt
  fi

done < missing_fnafaa_bacteria.txt

## These should be redownloaded
cat missing_fnafaa_bacteria_validation.txt | grep valid
```

We can manually download these with:

```bash
for i in $(cat missing_fnafaa_bacteria_validation.txt | grep valid | cut -d ' ' -f 2); do
  dirname=$(echo $i | rev | cut -d '/' -f 2 | rev | cut -d '_' -f 1-2 )
  mkdir -p missing/assembly_bacteria/ncbi_dataset/data/$dirname
  wget $i -P missing/assembly_bacteria/ncbi_dataset/data/$dirname/
done

cd
find -name '*.gz' -type f -exec gunzip {} \;
```

Retrieve the lists of re-downloaded files and place them in the correct directory.

```bash
cd missing/
find -name '*.gz' -type f -exec gunzip {} \;
find ~+ -type f -name "*.fna" > ../allfna_redownloaded.txt
find ~+ -type f -name "*.faa" > ../allfaa_redownloaded.txt
mv assembly_bacteria/ncbi_dataset/data/* ../assembly_bacteria/ncbi_dataset/data/
sed -i 's#missing/##g' ../allfna_redownloaded.txt
sed -i 's#missing/##g' ../allfaa_redownloaded.txt
cd ../
```

To reappend the re-downloaded files to the samplesheet, we can use the following script.

```bash
for type in  bacteria ; do
  while read line; do
      genomename=$(echo $line | cut -d ' ' -f 1)
      accnum=$(echo $line | cut -d ' ' -f 2)
      taxid=$(grep $accnum samplesheet_$type.txt | cut -f 2)
      fna=$(grep -e "$accnum" allfna_redownloaded.txt)
      faa=$(grep -e "$accnum" allfaa_redownloaded.txt)
      if [[ -z "$fna" && -z "$faa" ]]; then
        echo "missing BOTH: $genomename $accnum";
        echo "$genomename $accnum" >> missing_fnafaa_"$type"_redownloaded.txt
      elif [[ -z "$fna" && -n "$faa" ]]; then
        echo "missing FNA: $genomename $accnum"
        echo "$genomename $accnum" >> missing_fna_"$type"_redownloaded.txt
      elif [[ -n "$fna" && -z "$faa" ]]; then
        echo "missing FAA: $genomename $accnum"
        echo "$genomename $accnum" >> missing_faa_"$type"_redownloaded.txt
      else
        echo "${genomename}_${accnum},$taxid,$fna,$faa" >> samplesheet_"$type"_localpaths.csv
      fi
  done < <(cat missing_fnafaa_bacteria.txt)
done
```

Finally, we can combine the three samplesheets into one.

```bash
cat samplesheet_bacteria_localpaths.csv > samplesheet.csv
tail -n +2 samplesheet_viral_localpaths.csv >> samplesheet.csv
tail -n +2 samplesheet_archaea_localpaths.csv >> samplesheet.csv
```

TODO: DO SUBSAMPLING OF INITIAL CSV

-->
