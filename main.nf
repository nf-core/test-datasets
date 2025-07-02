#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GENERATE_EXAMPLE_GENOTYPES_VCFS } from './modules/generate_example_genotypes_vcfs.nf'
include { CHUNK_VCFS } from './modules/chunk_vcfs.nf'
include { CONCAT_CHUNKED_VCFS } from './modules/concat_chunked_vcfs.nf'
include { EXTRACT_SAMPLE_IDS } from './modules/extract_sample_ids.nf'
include { GENERATE_PHENO_COV } from './modules/generate_pheno_cov.nf'
include { INDEX_CHUNKED_VCFS } from './modules/index_chunked_vcfs.nf'
workflow {
    // Run the download process
    GENERATE_EXAMPLE_GENOTYPES_VCFS()

    def vcfs_with_chr = GENERATE_EXAMPLE_GENOTYPES_VCFS.out.vcfs
    	.flatten()
    	.map { file ->
        	def chr = file.name.toString().split("\\.")[1] // safer than `tokenize`
        	tuple(chr, file)
    	}

    // Feed the tuples into the chunking process
    CHUNK_VCFS(vcfs_with_chr)
    chr1_chunked = CHUNK_VCFS.out.chunked_vcfs.filter { chr, file -> chr == 'chr1'}.map { chr, file -> file }
    EXTRACT_SAMPLE_IDS(chr1_chunked)
    GENERATE_PHENO_COV(EXTRACT_SAMPLE_IDS.out.sample_ids)
    INDEX_CHUNKED_VCFS(CHUNK_VCFS.out.chunked_vcfs)
    all_chunked_vcfs = CHUNK_VCFS.out.chunked_vcfs.map { chr, file -> file }.collect()
    CONCAT_CHUNKED_VCFS(all_chunked_vcfs)
}
