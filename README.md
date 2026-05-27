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

## Files

- `testdata/faecalibacterium_taiwanense.fasta.gz` - _Faecalibacterium taiwanense_ genome assembly (gzipped FASTA)
- `samplesheet.tsv` - Input samplesheet for nf-core/bacmodel pipeline

## Usage

Use the test profile in nf-core/bacmodel:

```bash
nextflow run nf-core/bacmodel -profile test,docker
```

The test configuration automatically references this samplesheet from the test-datasets repository.

## Samplesheet Format

```tsv
sample	fasta
faecalibacterium_taiwanense	https://raw.githubusercontent.com/nf-core/test-datasets/bacmodel/testdata/faecalibacterium_taiwanense.fasta.gz
```
