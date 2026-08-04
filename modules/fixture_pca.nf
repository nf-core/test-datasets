// Principal components of the merged fixture genotypes, used as covariates.
// Exact (non-approximate) PCA, so the result is a deterministic function of the
// genotype matrix.
process FIXTURE_PCA {
	container "community.wave.seqera.io/library/plink2:2.0.0a.6.9--e6710830a4b7f0c6"

	input:
	tuple val(prefix), path(pgen), path(pvar), path(psam)

	output:
	path "${prefix}.eigenvec", emit: eigenvec
	path "${prefix}.eigenval", emit: eigenval

	script:
	"""
	plink2 \\
		--pfile ${prefix} \\
		--pca ${params.fixture_n_pcs_computed} \\
		--threads 1 \\
		--seed 1 \\
		--out ${prefix}
	"""

	stub:
	"""
	printf '#FID\\tIID\\tPC1\\n' > ${prefix}.eigenvec
	printf '0.0\\n' > ${prefix}.eigenval
	"""
}
