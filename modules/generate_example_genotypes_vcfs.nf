process GENERATE_EXAMPLE_GENOTYPES_VCFS {
        container "community.wave.seqera.io/library/bcftools_tabix_pip_tools:48085064a9189d8c"
        publishDir params.outdir_vcfs, mode: 'symlink'

        output:
        path "ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", emit: chr1_vcf
        path "ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", emit: chr2_vcf
        path "ALL.chr3.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", emit: chr3_vcf
        path "ALL.chr4.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", emit: chr4_vcf
        path "ALL.chr5.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", emit: chr5_vcf
        path "ALL.chr6.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", emit: chr6_vcf
        path "ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", emit: chr7_vcf
        path "ALL.chr8.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", emit: chr8_vcf
        path "ALL.chr9.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", emit: chr9_vcf
        path "ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", emit: chr10_vcf
        path "ALL.chr11.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", emit: chr11_vcf
        path "ALL.chr12.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", emit: chr12_vcf
        path "ALL.chr13.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", emit: chr13_vcf
        path "ALL.chr14.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", emit: chr14_vcf
        path "ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", emit: chr15_vcf
        path "ALL.chr16.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", emit: chr16_vcf
        path "ALL.chr17.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", emit: chr17_vcf
        path "ALL.chr18.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", emit: chr18_vcf
        path "ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", emit: chr19_vcf
        path "ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", emit: chr20_vcf
        path "ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", emit: chr21_vcf
        path "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", emit: chr22_vcf
        path "ALL.chrX.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", emit: chrX_vcf

	script:
	"""
	wget 
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	wget 
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	wget 
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	wget 
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	wget 
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	wget 
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	wget 
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	wget 
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	wget 
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	wget 
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	wget 
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	wget 
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	wget 
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	wget 
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	wget 
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	wget 
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	wget 
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	wget 
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	wget 
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	wget 
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	wget 
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	wget 
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	wget 
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	"""
