process INDEX_CHUNKED_VCFS {
	container "community.wave.seqera.io/library/bcftools_tabix_pip_tools:48085064a9189d8c"
	publishDir params.outdir_chunked_vcfs, mode: 'copy'

	input:
	tuple val(chr), path(vcf_file)

	output:
	tuple val(chr), path("*.vcf.gz.tbi"), emit: indexed_vcfs

	script:
	"""
	tabix -p vcf ${vcf_file}
	"""
}
