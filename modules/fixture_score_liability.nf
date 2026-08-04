// Compute each sample's genetic value for both simulated traits directly from the
// fixture genotypes, so the simulated signal is carried by the real genotype matrix
// rather than added as noise.
process FIXTURE_SCORE_LIABILITY {
	container "community.wave.seqera.io/library/plink2:2.0.0a.6.9--e6710830a4b7f0c6"

	input:
	tuple val(prefix), path(pgen), path(pvar), path(psam)
	path weights

	output:
	path "${prefix}.sscore", emit: sscore

	script:
	"""
	plink2 \\
		--pfile ${prefix} \\
		--score ${weights} 1 2 header-read cols=+scoresums \\
		--score-col-nums 3-4 \\
		--threads 1 \\
		--seed 1 \\
		--out ${prefix}
	"""

	stub:
	"""
	printf '#FID\\tIID\\tALLELE_CT\\tNAMED_ALLELE_DOSAGE_SUM\\tBETA_QT_AVG\\tBETA_BT_AVG\\tBETA_QT_SUM\\tBETA_BT_SUM\\n' > ${prefix}.sscore
	"""
}
