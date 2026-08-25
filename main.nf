#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GENERATE_GWAS_FIXTURES } from './modules/generate_gwas_fixtures.nf'

workflow {
    GENERATE_GWAS_FIXTURES()
}
