# Broken Samplesheets for input schema testing

The following input FASTQ samplesheets files are all broken in ways that are described in the file name.
Each one should give a valid and understandable error message as reported by nf-schema.
The first file `01` corresponds to a valid working sheet, thus not included here.

- `02-samplesheet_v3_doublesinglestrandednotuniquesampleid.tsv`
- `03-samplesheet_v3_wrongheadername.tsv`
- `04-samplesheet_v3_missingcolumns.tsv`
- `05-samplesheet_v3_sampleidmissing.tsv`
- `06-samplesheet_v3_libidmissing.tsv`
- `07-samplesheet_v3_lanenotnumeric.tsv`
- `08-samplesheet_v3_wrongcolourchem.tsv`
- `09-samplesheet_v3_wrongpairment.tsv`
- `10-samplesheet_v3_wrongstrandedness.tsv`
- `11-samplesheet_v3_wrongdamagetreatment.tsv`
- `12-samplesheet_v3_pathwithspaces.tsv`
- `13-samplesheet_v3_invalidfilextensions.tsv`
- `14-samplesheet_v3_bamnotsingleend.tsv`
- `15-samplesheet_v3_bammissingrefid.tsv`
- `16-samplesheet_v3_fastqwithrefid.tsv`
- `17-samplesheet_v3_r2setassingleend.tsv`
- `18-samplesheet_v3_r2missingr1.tsv`
- `19-samplesheet_v3_duplicatelibidinsample.tsv`
- `20-samplesheet_v3_pairedfastqadnbamsameline.tsv`
- `21-samplesheet_v3_singlefastqadnbamsameline.tsv`
- `22-samplesheet_v3_missingfile.tsv`
- `23-samplesheet_v3_missingr1orbam.tsv`
