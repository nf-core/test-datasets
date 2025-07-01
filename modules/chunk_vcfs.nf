process CHUNK_VCFS {
	container "community.wave.seqera.io/library/bcftools_tabix_pip_tools:48085064a9189d8c"
	publishDir params.outdir_chunked_vcfs, mode: 'copy'

	input:
	tuple val(chr), path(vcfs)

	output:
	path("${chr}_chunked.vcf.gz"), emit: chunked_vcfs

	script:
	"""
    # Proper VCF chunking with bcftools
    bcftools view -H ${vcfs} | head -n 4500 > variants.txt
    bcftools view -h ${vcfs} > header.txt
    cat header.txt variants.txt | bgzip > ${chr}_chunked.vcf.gz
	"""
}
