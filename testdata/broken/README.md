# Broken input sheets schema testing

## Samplesheets

The following input FASTQ samplesheets files are all broken in ways that are described in the file name.
Each one should give a valid and understandable error message as reported by nf-schema.
The first file `01` corresponds to a valid working sheet, thus not included here.

- `samplesheet/02-samplesheet_v3_doublesinglestrandednotuniquesampleid.tsv`
- `samplesheet/03-samplesheet_v3_wrongheadername.tsv`
- `samplesheet/04-samplesheet_v3_missingcolumns.tsv`
- `samplesheet/05-samplesheet_v3_sampleidmissing.tsv`
- `samplesheet/06-samplesheet_v3_libidmissing.tsv`
- `samplesheet/07-samplesheet_v3_lanenotnumeric.tsv`
- `samplesheet/08-samplesheet_v3_wrongcolourchem.tsv`
- `samplesheet/09-samplesheet_v3_wrongpairment.tsv`
- `samplesheet/10-samplesheet_v3_wrongstrandedness.tsv`
- `samplesheet/11-samplesheet_v3_wrongdamagetreatment.tsv`
- `samplesheet/12-samplesheet_v3_pathwithspaces.tsv`
- `samplesheet/13-samplesheet_v3_invalidfilextensions.tsv`
- `samplesheet/14-samplesheet_v3_bamnotsingleend.tsv`
- `samplesheet/15-samplesheet_v3_bammissingrefid.tsv`
- `samplesheet/16-samplesheet_v3_fastqwithrefid.tsv`
- `samplesheet/17-samplesheet_v3_r2setassingleend.tsv`
- `samplesheet/18-samplesheet_v3_r2missingr1.tsv`
- `samplesheet/19-samplesheet_v3_duplicatelibidinsample.tsv`
- `samplesheet/20-samplesheet_v3_pairedfastqadnbamsameline.tsv`
- `samplesheet/21-samplesheet_v3_singlefastqadnbamsameline.tsv`
- `samplesheet/22-samplesheet_v3_missingfile.tsv`
- `samplesheet/23-samplesheet_v3_missingr1orbam.tsv`

## Reference Sheets

The following input reference sheets files are all broken in ways that are described in the file name.
Each one should give a valid and understandable error message as reported by nf-schema.
The first file `01` corresponds to a valid working sheet, thus not included here.

- `referencesheet/01-reference_sheet_multiref_valid.csv`
- `referencesheet/02-reference_sheet_multiref_missingreferencename.csv`
- `referencesheet/03-reference_sheet_multiref_invalidreferencenamespaces.csv`
- `referencesheet/04-reference_sheet_multiref_missingfasta.csv`
- `referencesheet/05-reference_sheet_multiref_wrongfastasuffix.csv`
- `referencesheet/06-reference_sheet_multiref-nonuniquefasta.csv`
- `referencesheet/07-reference_sheet_multiref-indexfiletypoinpath.csv`
- `referencesheet/08-reference_sheet_multiref-badformattedcirculartarget.csv`
- `referencesheet/09-reference_sheet_multiref-pileupcallermissingdependency.csv`