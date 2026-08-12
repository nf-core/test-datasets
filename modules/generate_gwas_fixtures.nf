process GENERATE_GWAS_FIXTURES {
	container "community.wave.seqera.io/library/bcftools_tabix_pip_tools:48085064a9189d8c"
	publishDir "${params.outdir_fixtures}/genotypes", mode: 'copy', overwrite: true, pattern: "${params.fixture_prefix}_all.vcf.gz"
	publishDir "${params.outdir_fixtures}/pheno_cov", mode: 'copy', overwrite: true, pattern: "${params.fixture_prefix}.{pheno,qcovar,catcovar}"

	output:
	path "${params.fixture_prefix}_all.vcf.gz", emit: vcf
	path "${params.fixture_prefix}.pheno", emit: pheno
	path "${params.fixture_prefix}.qcovar", emit: qcovar
	path "${params.fixture_prefix}.catcovar", emit: catcovar

	script:
	"""
	python3 ${projectDir}/bin/generate_gwas_fixtures.py \
		--prefix ${params.fixture_prefix} \
		--chromosomes ${params.fixture_chromosomes} \
		--samples ${params.fixture_n_samples} \
		--variants-per-chromosome ${params.fixture_variants_per_chromosome} \
		--cases ${params.fixture_n_cases} \
		--seed ${params.fixture_seed}

	bgzip --threads 1 ${params.fixture_prefix}_all.vcf

	python3 ${projectDir}/bin/validate_gwas_fixtures.py \
		--vcf ${params.fixture_prefix}_all.vcf.gz \
		--pheno ${params.fixture_prefix}.pheno \
		--qcovar ${params.fixture_prefix}.qcovar \
		--catcovar ${params.fixture_prefix}.catcovar \
		--samples ${params.fixture_n_samples} \
		--chromosomes ${params.fixture_chromosomes} \
		--variants-per-chromosome ${params.fixture_variants_per_chromosome} \
		--cases ${params.fixture_n_cases}
	"""
}
