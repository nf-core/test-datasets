#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GENERATE_EXAMPLE_GENOTYPES_VCFS } from './modules/generate_example_genotypes_vcfs.nf'
include { CHUNK_VCFS } from './modules/chunk_vcfs.nf'
include { CONCAT_CHUNKED_VCFS } from './modules/concat_chunked_vcfs.nf'
include { EXTRACT_SAMPLE_IDS } from './modules/extract_sample_ids.nf'
include { GENERATE_PHENO_COV } from './modules/generate_pheno_cov.nf'
include { INDEX_CHUNKED_VCFS } from './modules/index_chunked_vcfs.nf'
include { FIXTURE_SELECT_SAMPLES } from './modules/fixture_select_samples.nf'
include { FIXTURE_THIN_CHROMOSOME } from './modules/fixture_thin_chromosome.nf'
include { FIXTURE_MERGE_CHROMOSOMES } from './modules/fixture_merge_chromosomes.nf'
include { FIXTURE_PCA } from './modules/fixture_pca.nf'
include { FIXTURE_CAUSAL_WEIGHTS } from './modules/fixture_causal_weights.nf'
include { FIXTURE_SCORE_LIABILITY } from './modules/fixture_score_liability.nf'
include { FIXTURE_PHENO_COVAR } from './modules/fixture_pheno_covar.nf'
include { FIXTURE_EXPORT_FORMATS } from './modules/fixture_export_formats.nf'
include { FIXTURE_COMPRESS_VCF } from './modules/fixture_compress_vcf.nf'
workflow {
    // Source-data stage: download the 1000 Genomes release and chunk it. Skipped
    // with --skip_source_generation, because its output is already committed under
    // results/chunked_vcfs and re-running it re-downloads tens of gigabytes.
    if (!params.skip_source_generation) {
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

    // GWAS fixture stage: derive the pipeline-level test fixtures from the
    // already-committed chromosome-chunked VCFs. Nothing here reads the network.
    def fixture_chrs = params.fixture_chromosomes.tokenize(',').collect { chr -> chr.trim() }
    def fixture_sources = channel
        .fromList(fixture_chrs)
        .map { chr -> tuple(chr, file("${params.chunked_vcfs_dir}/${chr}_chunked.vcf.gz", checkIfExists: true)) }

    // The sample subset is taken from the first fixture chromosome; every chunked
    // VCF carries the same 2504 samples in the same order.
    FIXTURE_SELECT_SAMPLES(
        file("${params.chunked_vcfs_dir}/${fixture_chrs.first()}_chunked.vcf.gz", checkIfExists: true)
    )

    FIXTURE_THIN_CHROMOSOME(
        fixture_sources,
        FIXTURE_SELECT_SAMPLES.out.samples,
        FIXTURE_SELECT_SAMPLES.out.sex
    )

    FIXTURE_MERGE_CHROMOSOMES(
        fixture_chrs,
        FIXTURE_THIN_CHROMOSOME.out.pgen
            .map { prefix, pgen, pvar, psam -> [pgen, pvar, psam] }
            .collect()
    )

    FIXTURE_PCA(FIXTURE_MERGE_CHROMOSOMES.out.pgen)

    FIXTURE_CAUSAL_WEIGHTS(
        FIXTURE_MERGE_CHROMOSOMES.out.pgen.map { prefix, pgen, pvar, psam -> pvar }
    )

    FIXTURE_SCORE_LIABILITY(
        FIXTURE_MERGE_CHROMOSOMES.out.pgen,
        FIXTURE_CAUSAL_WEIGHTS.out.weights
    )

    FIXTURE_PHENO_COVAR(
        FIXTURE_SCORE_LIABILITY.out.sscore,
        FIXTURE_PCA.out.eigenvec,
        FIXTURE_SELECT_SAMPLES.out.sex
    )

    // Per-chromosome and combined bundles are both emitted: leave-one-chromosome-out
    // methods need the per-chromosome split, single-cohort methods want the merge.
    def fixture_bundles = FIXTURE_THIN_CHROMOSOME.out.pgen.mix(FIXTURE_MERGE_CHROMOSOMES.out.pgen)

    FIXTURE_EXPORT_FORMATS(fixture_bundles)
    FIXTURE_COMPRESS_VCF(FIXTURE_EXPORT_FORMATS.out.vcf)
}
