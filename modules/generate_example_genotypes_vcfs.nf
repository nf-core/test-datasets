process GENERATE_EXAMPLE_GENOTYPES_VCFS {
        container "community.wave.seqera.io/library/bcftools_tabix_pip_tools:48085064a9189d8c"
        publishDir params.outdir_vcfs, mode: 'symlink'

        output:
        path "*.vcf.gz", emit: vcfs

	script:
	"""
    	for chr in {1..22}; do
        	fname="ALL.chr\${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
       		curl -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/\$fname
    	done
	curl -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz
	"""
}
