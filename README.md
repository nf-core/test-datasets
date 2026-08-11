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

## rnastructurome pipeline test data

Test data for the [nf-core/rnastructurome](https://github.com/nf-core/rnastructurome) pipeline lives on the `rnastructurome` branch under `testdata/`.

There are three independent test datasets, each exercising a different route through the pipeline:

- **`test_genome`**: tests the genome-alignment route. Run with `-profile test,docker`.
- **`test_transcriptome`**: tests the transcriptome-alignment route, including the `structextract` module. Run with `-profile test_transcriptome,docker`.
- **`test_prokaryote`**: tests the `jackknife` and `eval` modules using prokaryote 16S rRNA data and a known reference structure for that rRNA. Run with `-profile test_prokaryote,conda`.

**Human (HEK293) data** — the `test_genome` and `test_transcriptome` datasets are both derived from a real HEK293 dataset (accessions listed below). Reads aligning to the mitochondrial gene RNR1 (`ENST00000389680.2`) were extracted from the full dataset to create this small test data.

- For the genome route, the full GRCh38 FASTA and GTF were each subsetted to just the MT chromosome.
- For the transcriptome route, the RNR1 transcript was extracted from the full GRCh38 transcriptome FASTA, and the same MT GTF was reused.

**Prokaryote (E. coli) data** — reads aligning to 16S rRNA were extracted from a real _E. coli_ dataset (accession listed below) to create the tiny prokaryote dataset. The reference structure for 16S rRNA was supplied by Danny Incarnato and is used by both the `jackknife` and `eval` modules. The FASTA and GTF files were created from the 16S rRNA reference structure.

### Directory structure

```
testdata/
├── HEK293T_untreated_r1.fastq.gz       # shared SHAPE-seq reads (untreated, GSM4333255)
├── HEK293T_treated_r1.fastq.gz         # shared SHAPE-seq reads (treated, GSM4333256)
├── test_genome/
│   ├── Homo_sapiens.GRCh38.MT.fa       # human mitochondrial chromosome
│   ├── Homo_sapiens.GRCh38.MT.gtf      # MT genome annotation
│   └── samplesheet.genome_test.csv
├── test_transcriptome/
│   ├── Homo_sapiens.GRCh38.ENST00000389680.fa  # single-transcript FASTA (MT-RNR1)
│   ├── Homo_sapiens.GRCh38.MT.gtf              # MT genome annotation
│   └── samplesheet.transcriptome_test.csv
└── test_prokaryote/
    ├── 16S_rRNA.fa                      # E. coli 16S rRNA FASTA
    ├── 16S_rRNA.gtf                     # matching single-transcript GTF
    ├── 16S_rRNA.reference.db            # dot-bracket reference structure for rf-eval/rf-jackknife
    ├── E_coli_DH5a_treated.16S_rRNA.fastq.gz  # treated, GSM7885842
    └── samplesheet.prokaryote_test.csv
```

## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).

[^1]: From [stackoverflow](https://stackoverflow.com/a/60846265/11502856)
