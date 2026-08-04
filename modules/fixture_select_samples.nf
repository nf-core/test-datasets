// Deterministically select the sample subset used by every GWAS fixture file.
//
// Selection is a fixed stride over the lexicographically sorted sample IDs of the
// source VCF, so the same subset is produced on every run without any RNG.
// The source VCFs carry no sex information (the committed chrX chunk lies entirely
// within PAR1, where both sexes are diploid), so sex is assigned deterministically
// by alternating over the sorted selection.
process FIXTURE_SELECT_SAMPLES {
	container "community.wave.seqera.io/library/bcftools_tabix_pip_tools:48085064a9189d8c"
	publishDir params.outdir_fixture_pheno_cov, mode: 'copy'

	input:
	path vcf_file

	output:
	path "fixture_samples.tsv", emit: samples
	path "fixture_sex.tsv", emit: sex

	script:
	"""
	bcftools query -l ${vcf_file} \\
		| LC_ALL=C sort \\
		| awk 'NR % ${params.fixture_sample_stride} == 1' \\
		| head -n ${params.fixture_n_samples} > selected_ids.txt

	selected=\$(wc -l < selected_ids.txt)
	if [ "\${selected}" -ne "${params.fixture_n_samples}" ]; then
		echo "ERROR: selected \${selected} samples, expected ${params.fixture_n_samples}" >&2
		exit 1
	fi

	printf '#FID\\tIID\\n' > fixture_samples.tsv
	awk 'BEGIN { OFS = "\\t" } { print \$1, \$1 }' selected_ids.txt >> fixture_samples.tsv

	printf '#FID\\tIID\\tSEX\\n' > fixture_sex.tsv
	awk 'BEGIN { OFS = "\\t" } { print \$1, \$1, (NR % 2 == 1 ? 1 : 2) }' selected_ids.txt >> fixture_sex.tsv
	"""

	stub:
	"""
	printf '#FID\\tIID\\n' > fixture_samples.tsv
	printf '#FID\\tIID\\tSEX\\n' > fixture_sex.tsv
	for i in 1 2 3; do
		printf 'SAMPLE%s\\tSAMPLE%s\\n' "\${i}" "\${i}" >> fixture_samples.tsv
		printf 'SAMPLE%s\\tSAMPLE%s\\t1\\n' "\${i}" "\${i}" >> fixture_sex.tsv
	done
	"""
}
