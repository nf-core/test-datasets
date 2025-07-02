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
    # Create deterministic VCF files for testing
    for chr in 1 X; do
        if [ "\$chr" = "X" ]; then
            fname="ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz"
        else
            fname="ALL.chr\${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
        fi
        
        # Create VCF with proper 1000 Genomes-like header
        cat > \${fname%.gz} << EOF
    ##fileformat=VCFv4.1
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##contig=<ID=chr\${chr}>
    ##reference=file:///path/to/human_g1k_v37.fasta
    ##source=1000GenomesPhase3Pipeline
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096	HG00097	HG00099
    EOF

        # Add deterministic variants (same every time)
        for i in \$(seq 1 4600); do
            pos=\$((i * 10000))
            
            # Use deterministic patterns instead of random
            case \$((i % 4)) in
                0) ref="A"; alt="T" ;;
                1) ref="T"; alt="G" ;;
                2) ref="G"; alt="C" ;;
                3) ref="C"; alt="A" ;;
            esac
            
            # Deterministic genotypes based on position
            case \$((i % 3)) in
                0) gt1="0|0"; gt2="0|1"; gt3="1|1" ;;
                1) gt1="0|1"; gt2="1|0"; gt3="0|0" ;;
                2) gt1="1|1"; gt2="0|0"; gt3="0|1" ;;
            esac
            
            echo -e "chr\${chr}\t\${pos}\trs\${i}\t\${ref}\t\${alt}\t100\tPASS\tAC=2;AF=0.5;AN=6;DP=50\tGT\t\${gt1}\t\${gt2}\t\${gt3}" >> \${fname%.gz}
        done
        
        # Compress the VCF
        bgzip \${fname%.gz}
        echo "Created \$fname"
        
        echo "Generated deterministic \$fname with 4600 variants"
    done
    
    # Copy or symlink the files to params.outdir_vcfs
    mkdir -p ${params.outdir_vcfs}
    cp ALL.chr*.vcf.gz ${params.outdir_vcfs}/
    """
}