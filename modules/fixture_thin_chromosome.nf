// Build the per-chromosome PLINK 2 fixture from an already-committed chunked VCF.
//
// The sample subset is fixed by FIXTURE_SELECT_SAMPLES, variants are restricted to
// biallelic ACGT SNPs above the fixture MAF floor, and the survivors are LD-thinned
// with plink2 --indep-pairwise. No step here reads the network.
process FIXTURE_THIN_CHROMOSOME {
	container "community.wave.seqera.io/library/plink2:2.0.0a.6.9--e6710830a4b7f0c6"
	publishDir params.outdir_fixture_genotypes, mode: 'copy', pattern: "*.{pgen,pvar,psam}"

	input:
	tuple val(chr), path(vcf_file)
	path samples
	path sex

	output:
	tuple val("${params.fixture_prefix}_${chr}"), path("${params.fixture_prefix}_${chr}.pgen"), path("${params.fixture_prefix}_${chr}.pvar"), path("${params.fixture_prefix}_${chr}.psam"), emit: pgen

	script:
	def prefix = "${params.fixture_prefix}_${chr}"
	"""
	# Sample subset plus basic variant QC. INFO/QUAL/FILTER are dropped from the
	# .pvar (pvar-cols=xheader): the source INFO fields describe all 2504 samples
	# and would be stale for the 500-sample fixture.
	plink2 \\
		--vcf ${vcf_file} \\
		--double-id \\
		--keep ${samples} \\
		--update-sex ${sex} \\
		--snps-only just-acgt \\
		--max-alleles 2 \\
		--set-missing-var-ids '@:#:\$r:\$a' \\
		--rm-dup exclude-all \\
		--maf ${params.fixture_maf} \\
		--threads 1 \\
		--seed 1 \\
		--make-pgen pvar-cols=xheader \\
		--out filtered

	plink2 \\
		--pfile filtered \\
		--indep-pairwise ${params.fixture_ld_window} ${params.fixture_ld_step} ${params.fixture_ld_r2} \\
		--threads 1 \\
		--seed 1 \\
		--out pruned

	plink2 \\
		--pfile filtered \\
		--extract pruned.prune.in \\
		--threads 1 \\
		--seed 1 \\
		--make-pgen pvar-cols=xheader \\
		--out ${prefix}
	"""

	stub:
	def prefix = "${params.fixture_prefix}_${chr}"
	"""
	touch ${prefix}.pgen
	printf '#CHROM\\tPOS\\tID\\tREF\\tALT\\n' > ${prefix}.pvar
	printf '#FID\\tIID\\tSEX\\n' > ${prefix}.psam
	"""
}
