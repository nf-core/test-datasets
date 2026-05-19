# nf-core/bacmodel Test Dataset

## Overview

This directory contains test data for the nf-core/bacmodel pipeline.

## Test Genome

**Species:** _Faecalibacterium taiwanense_  
**Source:** BIP-2026-46 Silva Assemblies (IKMB, sample 3606MB_1118_055_H10)  
**Assembly:** Bacterial isolate genome assembly  
**Size:** ~2.7 Mbp (13 contigs)  
**File:** `faecalibacterium_taiwanense.fasta.gz`

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
faecalibacterium_taiwanense	https://raw.githubusercontent.com/nf-core/test-datasets/bacmodel/bacmodel/testdata/faecalibacterium_taiwanense.fasta.gz
```
