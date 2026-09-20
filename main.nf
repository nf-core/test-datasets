#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GENERATE_GWAS_FIXTURES }   from './modules/generate_gwas_fixtures.nf'
include { PLINK2_GWAS_DERIVATIVES }  from './modules/plink2_gwas_derivatives.nf'
include { VALIDATE_GWAS_FIXTURES }   from './modules/validate_gwas_fixtures.nf'

workflow {
    GENERATE_GWAS_FIXTURES()
    PLINK2_GWAS_DERIVATIVES(GENERATE_GWAS_FIXTURES.out.vcf)
    VALIDATE_GWAS_FIXTURES(
        GENERATE_GWAS_FIXTURES.out.vcf,
        GENERATE_GWAS_FIXTURES.out.pheno,
        GENERATE_GWAS_FIXTURES.out.qcovar,
        GENERATE_GWAS_FIXTURES.out.catcovar,
        GENERATE_GWAS_FIXTURES.out.relational,
        PLINK2_GWAS_DERIVATIVES.out.pgen_all,
        PLINK2_GWAS_DERIVATIVES.out.pgen_chr,
        PLINK2_GWAS_DERIVATIVES.out.bed_all,
    )
}
