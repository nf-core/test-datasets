process GENERATE_EXAMPLE_GENOTYPES_VCFS {
    container "community.wave.seqera.io/library/bcftools_tabix_pip_tools:48085064a9189d8c"
    publishDir params.outdir_vcfs, mode: 'symlink'

    output:
    path "ALL.chr*.vcf.gz", emit: vcfs    
    
    script:
    """
    for chr in {1..22}; do
        fname="ALL.chr\${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
        curl -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/\$fname
    done
    curl -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz
    """

    stub:
    """
    echo "Creating mock VCF files..."

    # Create header content
    header='##fileformat=VCFv4.1
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1'

    # Create chromosomes 1-22
    for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do
        fname="ALL.chr\${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
        
        {
            echo "\$header"
            seq 1 4501 | awk -v chr=\$chr '{printf "%s\\t%d\\trs%d\\tA\\tT\\t60\\tPASS\\tAC=1\\tGT\\t0/1\\n", chr, \$1*1000, \$1}'
        } | gzip > \$fname
        
        echo "Created \$fname"
    done

    # Create X chromosome
    {
        echo "\$header"
        seq 1 4501 | awk '{printf "X\\t%d\\trs%d\\tA\\tT\\t60\\tPASS\\tAC=1\\tGT\\t0/1\\n", \$1*1000, \$1}'
    } | gzip > ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz

    echo "Created ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz"

    # Verify all files were created
    echo "Final file count:"
    ls -1 ALL.chr*.vcf.gz | wc -l
    ls -la ALL.chr*.vcf.gz
    """
}
