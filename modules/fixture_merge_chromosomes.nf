// Merge the per-chromosome PLINK 2 fixtures into one multi-chromosome bundle.
//
// The merge list is built from the ordered chromosome list rather than from the
// staged file order, so the merged variant order is reproducible.
process FIXTURE_MERGE_CHROMOSOMES {
	container "community.wave.seqera.io/library/plink2:2.0.0a.6.9--e6710830a4b7f0c6"
	publishDir params.outdir_fixture_genotypes, mode: 'copy', pattern: "*.{pgen,pvar,psam}"

	input:
	val chrs
	path psets

	output:
	tuple val("${params.fixture_prefix}_all"), path("${params.fixture_prefix}_all.pgen"), path("${params.fixture_prefix}_all.pvar"), path("${params.fixture_prefix}_all.psam"), emit: pgen

	script:
	def prefix = "${params.fixture_prefix}_all"
	def merge_list = chrs.collect { chr -> "${params.fixture_prefix}_${chr}" }.join('\\n')
	"""
	printf '${merge_list}\\n' > merge_list.txt

	plink2 \\
		--pmerge-list merge_list.txt pfile \\
		--threads 1 \\
		--seed 1 \\
		--out ${prefix}
	"""

	stub:
	def prefix = "${params.fixture_prefix}_all"
	"""
	touch ${prefix}.pgen
	printf '#CHROM\\tPOS\\tID\\tREF\\tALT\\n' > ${prefix}.pvar
	printf '#FID\\tIID\\tSEX\\n' > ${prefix}.psam
	"""
}
