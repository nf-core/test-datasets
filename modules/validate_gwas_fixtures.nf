// Semantic validation of the complete published family, run after the PLINK derivatives exist so the
// cohort manifest's PLINK 2 selection can be resolved against real files rather than assumed.
process VALIDATE_GWAS_FIXTURES {
	container "community.wave.seqera.io/library/bcftools_tabix_pip_tools:48085064a9189d8c"

	input:
	path vcf
	path pheno
	path qcovar
	path catcovar
	path relational_dir
	tuple path(pgen), path(psam), path(pvar)
	tuple path(bed), path(bim), path(fam)

	output:
	path "fixture_summary.tsv", emit: summary

	script:
	"""
	python3 ${projectDir}/bin/validate_gwas_fixtures.py \
		--vcf ${vcf} \
		--pheno ${pheno} \
		--qcovar ${qcovar} \
		--catcovar ${catcovar} \
		--pgen ${pgen} \
		--psam ${psam} \
		--pvar ${pvar} \
		--bed ${bed} \
		--bim ${bim} \
		--fam ${fam} \
		--relational-dir ${relational_dir} \
		--samples ${params.fixture_n_samples} \
		--chromosomes ${params.fixture_chromosomes} \
		--variants-per-chromosome ${params.fixture_variants_per_chromosome} \
		--cases ${params.fixture_n_cases} \
		| tee fixture_summary.tsv
	"""
}
