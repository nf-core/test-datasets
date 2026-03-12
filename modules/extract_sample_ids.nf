process EXTRACT_SAMPLE_IDS {
	container "community.wave.seqera.io/library/r-base:4.4.3--1e564c44feffeaa0"
        publishDir params.outdir_pheno_cov, mode: 'symlink'

	input:
	path vcf_file

	output:
	path "sample_ids.txt", emit: sample_ids

	script:
	"""
    	zcat $vcf_file | grep '#CHROM' | cut -f10- | tr '\t' '\n' > sample_ids.txt
	"""
}
