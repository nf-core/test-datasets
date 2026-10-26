# nf-core/bacmodel Test Dataset

## Overview

This directory contains test data for the nf-core/bacmodel pipeline, which performs functional annotation and metabolic modeling of bacterial genomes.

The test dataset consists of a complete bacterial genome assembly from _Faecalibacterium taiwanense_, a beneficial gut bacterium known for its anti-inflammatory properties and butyrate production. This species is an excellent test case as it contains diverse functional features including secondary metabolite biosynthesis, carbohydrate-active enzymes, and well-characterized metabolic pathways - enabling comprehensive testing of all bacmodel annotation and modeling tools.

## Sampling Procedure

The test genome was selected from the Global Microbiome Conservancy (GMbC)

### Test Genome

**Species:** _Faecalibacterium taiwanense_  
**Source:** Global Microbiome Conservancy (GMbC, https://microbiomeconservancy.org/)
**Assembly:** Bacterial isolate genome assembly  
**Size:** ~2.7 Mbp (13 contigs)  
**File:** `faecalibacterium_taiwanense.fasta.gz`

### Assembly Details

**Sequencing Platform:** Illumina paired-end sequencing  
**Assembly Method:** Shovill (via Bactopia pipeline, https://github.com/bactopia/bactopia)

The genome assembly was obtained from the Global Microbiome Conservancy (GMbC) collection, representing a _Faecalibacterium taiwanense_ isolate from a Nigeria rural sample. This assembly was selected for testing due to its compact size (802 KB compressed, 13 contigs, ~2.7 Mbp) while maintaining complete functional features.

**No subsampling or downsampling was applied.** The complete assembly is used as-is to ensure all annotation and modeling tools have sufficient data for meaningful analysis while maintaining fast test execution times.

**Data Availability:** This genome is redistributed with permission as part of the nf-core test-datasets collection for reproducible pipeline testing.

## Multi-Sample Dataset

The single-genome dataset above cannot exercise any pipeline code path that only
differs between one sample and several (e.g. how per-sample outputs are staged
and collected into `BACMODEL_SUMMARY`). This second dataset adds three more
gut isolate genomes, chosen to span distinct genera/phyla rather than being
near-duplicate strains, so that multi-sample runs are also a meaningful test
of community-relevant diversity (useful groundwork for future community-level
modeling steps, e.g. MICOM/ecoGS, which need metabolically distinct members).

| Species | Genus | Phylum | Country | Size | Contigs | File |
|---|---|---|---|---|---|---|
| _Prevotella copri_A_ | Prevotella | Bacteroidota | Ghana | ~3.9 Mbp | 95 | `prevotella_copri.fasta.gz` |
| _Faecalibacterium duncaniae_ | Faecalibacterium | Bacillota_A | Ghana | ~2.8 Mbp | 60 | `faecalibacterium_duncaniae.fasta.gz` |
| _Dorea formicigenerans_ | Dorea | Bacillota_A | Tanzania | ~3.0 Mbp | 65 | `dorea_formicigenerans.fasta.gz` |

**Source:** Global Microbiome Conservancy (GMbC, https://microbiomeconservancy.org/)

**Samplesheet:** `samplesheet_multisample.tsv` lists these three samples only
(the original `faecalibacterium_taiwanense` sample is not included, to keep
this dataset's runtime independent of the single-sample one).

## Files

- `testdata/faecalibacterium_taiwanense.fasta.gz` - _Faecalibacterium taiwanense_ genome assembly (gzipped FASTA)
- `testdata/prevotella_copri.fasta.gz` - _Prevotella copri_A_ genome assembly (gzipped FASTA)
- `testdata/faecalibacterium_duncaniae.fasta.gz` - _Faecalibacterium duncaniae_ genome assembly (gzipped FASTA)
- `testdata/dorea_formicigenerans.fasta.gz` - _Dorea formicigenerans_ genome assembly (gzipped FASTA)
- `samplesheet.tsv` - Single-sample input samplesheet for nf-core/bacmodel pipeline
- `samplesheet_multisample.tsv` - Three-sample input samplesheet for nf-core/bacmodel pipeline

## Usage

Use the test profile in nf-core/bacmodel:

```bash
nextflow run nf-core/bacmodel -profile test,docker
```

For the multi-sample test profile:

```bash
nextflow run nf-core/bacmodel -profile test_multisample,docker
```

The test configurations automatically reference the relevant samplesheet from the test-datasets repository.

## Samplesheet Format

```tsv
sample	fasta
faecalibacterium_taiwanense	https://raw.githubusercontent.com/nf-core/test-datasets/bacmodel/testdata/faecalibacterium_taiwanense.fasta.gz
```

```tsv
sample	fasta	medium_gapseq	medium_gapseq_csv	medium_carveme	medium_carveme_tsv
prevotella_copri	https://raw.githubusercontent.com/nf-core/test-datasets/bacmodel/testdata/prevotella_copri.fasta.gz			LB
faecalibacterium_duncaniae	https://raw.githubusercontent.com/nf-core/test-datasets/bacmodel/testdata/faecalibacterium_duncaniae.fasta.gz			LB
dorea_formicigenerans	https://raw.githubusercontent.com/nf-core/test-datasets/bacmodel/testdata/dorea_formicigenerans.fasta.gz			LB
```
