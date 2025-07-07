process GENERATE_PHENO_COV {
	container "community.wave.seqera.io/library/r-base:4.4.3--1e564c44feffeaa0"
	publishDir params.outdir_pheno_cov, mode: 'copy'

	input:
	path sample_ids

	output:
	path "example.pheno", emit: pheno
	path "example.covar", emit: covar

	script:
	"""
	#!/usr/bin/env Rscript
	#make a phenotype
	#Here, a not too bad tutorial on different techniques on how to simulate data
	# https://aosmith.rbind.io/2018/08/29/getting-started-simulating-data/

	#Here I used the blog's proposed way of simulating data for a regression analysis
	#We will use the generated data slightly different, but hopefully good enough to
	# actually get some results
    	ids <- readLines("${sample_ids}")
    	n <- length(ids)
	set.seed(16)
	y = rnorm(n = n, mean = 0, sd = 1)
	x1 = runif(n = n, min = 1, max = 2)
	x2 = runif(n = n, min = 200, max = 300)

	# Write this to one phenodata file and one covardata file 
	# first column, unique ids, second column family ids, remaining columns are 
	# phenotyp or covariate columns (here individual IDs are family IDs)
	example.pheno <- data.frame(ids=ids, fam=ids, pheno=y)
	example.covar <- data.frame(ids=ids, fam=ids, cov1=x1, cov2=x2)
	# Write to tab-delimited files without headers or row names
        write.table(example.pheno, file = "example.pheno", sep = "\\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        write.table(example.covar, file = "example.covar", sep = "\\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
	"""
}
