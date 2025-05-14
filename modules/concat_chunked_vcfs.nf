process CONCAT_CHUNKED_VCFS {
	container "community.wave.seqera.io/library/bcftools_tabix_pip_tools:48085064a9189d8c"
	publishDir params.outdir_chunked_vcfs, mode: 'copy'

	input:
	path vcf_files

	output:
	path "combined_chunked.vcf.gz"
	path "combined_chunked.vcf.gz.tbi"

	script:
	"""
	echo "VCFs to concat:" > concat_debug.txt
	ls -lh ${vcf_files} >> concat_debug.txt

	bcftools concat -Oz -o combined_chunked.vcf.gz ${vcf_files.join(' ')}
	tabix -p vcf combined_chunked.vcf.gz	
	"""
}
