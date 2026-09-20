// PLINK derivatives of the canonical compact VCF.
//
// The three conversions below reproduce, flag for flag, the ones nf-core/gwas performs when it
// prepares this fixture family for its own test routes, so the published bundles are the exact bytes
// that pipeline consumes. The container pins PLINK v2.0.0-a.6.9, the same build the pipeline's
// `plink2` modules declare; changing it invalidates the committed byte-identical derivatives.
process PLINK2_GWAS_DERIVATIVES {
	container "community.wave.seqera.io/library/plink2:2.0.0a.6.9--e6710830a4b7f0c6"
	cpus 2
	memory 6.GB
	publishDir "${params.outdir_fixtures}/genotypes", mode: 'copy', overwrite: true, pattern: "*.{pgen,psam,pvar,bed,bim,fam}"

	input:
	path vcf

	output:
	tuple path("${params.fixture_prefix}_all.pgen"), path("${params.fixture_prefix}_all.psam"), path("${params.fixture_prefix}_all.pvar"), emit: pgen_all
	tuple path("${params.fixture_prefix}_chr${params.fixture_derivative_chromosome}.pgen"), path("${params.fixture_prefix}_chr${params.fixture_derivative_chromosome}.psam"), path("${params.fixture_prefix}_chr${params.fixture_derivative_chromosome}.pvar"), emit: pgen_chr
	tuple path("${params.fixture_prefix}_all.bed"), path("${params.fixture_prefix}_all.bim"), path("${params.fixture_prefix}_all.fam"), emit: bed_all

	script:
	def prefix = params.fixture_prefix
	def chromosome = params.fixture_derivative_chromosome
	def threads = task.cpus
	def mem_mb = task.memory.toMega()
	"""
	plink2 \
		--vcf "${vcf}" \
		--double-id \
		--threads "${threads}" \
		--memory "${mem_mb}" \
		--make-pgen \
		--out "${prefix}_all"

	plink2 \
		--vcf "${vcf}" \
		--double-id --chr ${chromosome} \
		--threads "${threads}" \
		--memory "${mem_mb}" \
		--make-pgen \
		--out "${prefix}_chr${chromosome}"

	plink2 \
		--pfile "${prefix}_all" \
		--threads "${threads}" \
		--memory "${mem_mb}" \
		--make-bed \
		--out "${prefix}_all" \
		--hard-call-threshold ${params.fixture_hard_call_threshold}
	"""
}
