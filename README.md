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

## Broken Samplesheets

To help improve schema checking, we've taking then main `test.csv`, and added a few variants which have various errors.

Each file _should_ fail and give an error message from nf-schema.

- `samplesheets/broken/test_duplicate_id_and_path.csv`: has cells in both `id` and `fasta_aa` duplicated when they should be unique
- `samplesheets/broken/test_duplicate_id_only.csv`: has cells in only `id` duplicated, when all cells should be unique
- `samplesheets/broken/test_missing_both_paths.csv`: has a row where both required `fasta_dna` and `fasta_aa` paths are missing
- `samplesheets/broken/test_missing_required_column.csv`: missing the required `taxid` column
- `samplesheets/broken/test_non_existent_file.csv`: has a path to a `fasta_dna` filepath that doesn't exist
