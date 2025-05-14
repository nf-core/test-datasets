process CHUNK_VCFS {
	container "community.wave.seqera.io/library/bcftools_tabix_pip_tools:48085064a9189d8c"
	publishDir params.outdir_chunked_vcfs, mode: 'copy'

	input:
	tuple val(chr), path(vcfs)

	output:
	path("${chr}_chunked.vcf.gz"), emit: chunked_vcfs

	script:
	"""
	bcftools view ${vcfs} | awk 'BEGIN {h=1; n=4500} /^#/ {print; next} {if (h <= n) {print; h++}}' | bgzip >${chr}_chunked.vcf.gz
	tabix -p vcf ${chr}_chunked.vcf.gz
	"""
}
