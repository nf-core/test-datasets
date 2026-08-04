// Block-compress and index the exported fixture VCFs.
//
// The ##fileDate header line plink2 writes carries the run date, which would make
// the committed fixture differ on every regeneration; it is dropped so the VCF is
// byte-reproducible.
process FIXTURE_COMPRESS_VCF {
	container "community.wave.seqera.io/library/bcftools_tabix_pip_tools:48085064a9189d8c"
	publishDir params.outdir_fixture_genotypes, mode: 'copy'

	input:
	tuple val(prefix), path(vcf_file)

	output:
	tuple val(prefix), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), emit: vcf

	script:
	"""
	grep -v '^##fileDate=' ${vcf_file} | bgzip -c > ${prefix}.vcf.gz
	tabix -p vcf ${prefix}.vcf.gz
	"""

	stub:
	"""
	printf '##fileformat=VCFv4.3\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\n' | bgzip -c > ${prefix}.vcf.gz
	touch ${prefix}.vcf.gz.tbi
	"""
}
