#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {GENERATE_EXAMPLE_GENOTYPES_VCFS } from 
'./modules/generate_example_genotypes_vcfs.nf'

workflow {
	GENERATE_EXAMPLE_GENOTYPES_VCFS()
}
