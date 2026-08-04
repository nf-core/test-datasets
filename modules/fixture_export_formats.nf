// Re-encode a PLINK 2 fixture bundle as PLINK 1 binary and as VCF, so all three
// encodings describe exactly the same samples and variants.
process FIXTURE_EXPORT_FORMATS {
	container "community.wave.seqera.io/library/plink2:2.0.0a.6.9--e6710830a4b7f0c6"
	publishDir params.outdir_fixture_genotypes, mode: 'copy', pattern: "*.{bed,bim,fam}"

	input:
	tuple val(prefix), path(pgen), path(pvar), path(psam)

	output:
	tuple val(prefix), path("${prefix}.bed"), path("${prefix}.bim"), path("${prefix}.fam"), emit: bed
	tuple val(prefix), path("${prefix}.vcf"), emit: vcf

	script:
	"""
	plink2 \\
		--pfile ${prefix} \\
		--threads 1 \\
		--seed 1 \\
		--make-bed \\
		--out ${prefix}

	plink2 \\
		--pfile ${prefix} \\
		--threads 1 \\
		--seed 1 \\
		--export vcf id-paste=iid \\
		--out ${prefix}
	"""

	stub:
	"""
	touch ${prefix}.bed ${prefix}.bim ${prefix}.fam
	printf '##fileformat=VCFv4.3\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\n' > ${prefix}.vcf
	"""
}
