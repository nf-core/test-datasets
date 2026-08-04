// Assemble the headered phenotype and covariate fixtures from the genotype-derived
// scores, the genotype-derived principal components and the deterministic sex
// assignment.
//
// The quantitative trait and the binary trait's liability are built as
//   sqrt(h2) * standardised genetic score
//   + sqrt(c2) * standardised covariate effect
//   + sqrt(1 - h2 - c2) * noise
// so both traits carry real genetic signal from the fixture genotypes and both
// respond to the supplied covariates. The binary trait is the liability
// dichotomised at the exact empirical quantile matching the target prevalence, so
// the case count is fixed rather than sampled.
process FIXTURE_PHENO_COVAR {
	container "community.wave.seqera.io/library/r-base:4.4.3--1e564c44feffeaa0"
	publishDir params.outdir_fixture_pheno_cov, mode: 'copy'

	input:
	path sscore
	path eigenvec
	path sex

	output:
	path "${params.fixture_prefix}.pheno", emit: pheno
	path "${params.fixture_prefix}.covar", emit: covar
	path "${params.fixture_prefix}.qcovar", emit: qcovar
	path "${params.fixture_prefix}.catcovar", emit: catcovar
	path "fixture_summary.txt", emit: summary

	script:
	"""
	#!/usr/bin/env Rscript
	set.seed(${params.fixture_seed},
			 kind = "Mersenne-Twister",
			 normal.kind = "Inversion",
			 sample.kind = "Rejection")

	read_headed <- function(path) {
		frame <- read.delim(path, header = TRUE, check.names = FALSE,
							comment.char = "", stringsAsFactors = FALSE)
		colnames(frame)[1] <- sub("^#", "", colnames(frame)[1])
		frame
	}

	standardise <- function(x) (x - mean(x)) / sd(x)

	scores <- read_headed("${sscore}")
	pcs <- read_headed("${eigenvec}")
	sexes <- read_headed("${sex}")

	ids <- scores[["IID"]]
	n <- length(ids)
	pcs <- pcs[match(ids, pcs[["IID"]]), , drop = FALSE]
	sexes <- sexes[match(ids, sexes[["IID"]]), , drop = FALSE]
	stopifnot(!anyNA(pcs[["IID"]]), !anyNA(sexes[["IID"]]))

	h2 <- ${params.fixture_h2}
	c2 <- ${params.fixture_covariate_variance}
	e2 <- 1 - h2 - c2
	stopifnot(e2 > 0)

	genetic_qt <- standardise(scores[["BETA_QT_SUM"]])
	genetic_bt <- standardise(scores[["BETA_BT_SUM"]])

	sex_code <- as.integer(sexes[["SEX"]])
	age <- pmin(pmax(round(rnorm(n, mean = 55, sd = 8)), 18), 90)
	covariate_effect <- standardise(0.6 * (sex_code == 2) + 0.03 * (age - 55))

	quantitative <- sqrt(h2) * genetic_qt +
		sqrt(c2) * covariate_effect +
		sqrt(e2) * rnorm(n)

	liability <- sqrt(h2) * genetic_bt +
		sqrt(c2) * covariate_effect +
		sqrt(e2) * rnorm(n)

	n_cases <- round(${params.fixture_prevalence} * n)
	threshold <- sort(liability, decreasing = TRUE)[n_cases]
	binary <- ifelse(liability >= threshold, 2L, 1L)

	pc_names <- paste0("PC", seq_len(${params.fixture_n_pcs}))
	pc_values <- lapply(pc_names, function(name) round(pcs[[name]], 6))
	names(pc_values) <- pc_names

	pheno <- data.frame(FID = ids, IID = ids,
						QT = round(quantitative, 6), BT = binary,
						stringsAsFactors = FALSE)
	covar <- data.frame(c(list(FID = ids, IID = ids,
							   SEX = sex_code, AGE = age), pc_values),
						stringsAsFactors = FALSE)
	qcovar <- data.frame(c(list(FID = ids, IID = ids, AGE = age), pc_values),
						 stringsAsFactors = FALSE)
	catcovar <- data.frame(FID = ids, IID = ids, SEX = sex_code,
						   stringsAsFactors = FALSE)

	write_tsv <- function(frame, path) {
		write.table(frame, file = path, sep = "\\t", quote = FALSE,
					row.names = FALSE, col.names = TRUE)
	}
	write_tsv(pheno, "${params.fixture_prefix}.pheno")
	write_tsv(covar, "${params.fixture_prefix}.covar")
	write_tsv(qcovar, "${params.fixture_prefix}.qcovar")
	write_tsv(catcovar, "${params.fixture_prefix}.catcovar")

	report <- c(
		sprintf("samples\\t%d", n),
		sprintf("target_h2\\t%.3f", h2),
		sprintf("target_covariate_variance\\t%.3f", c2),
		sprintf("target_prevalence\\t%.3f", ${params.fixture_prevalence}),
		sprintf("cases\\t%d", sum(binary == 2L)),
		sprintf("controls\\t%d", sum(binary == 1L)),
		sprintf("cor_QT_genetic_score\\t%.4f", cor(quantitative, genetic_qt)),
		sprintf("cor_QT_covariate_effect\\t%.4f", cor(quantitative, covariate_effect)),
		sprintf("cor_liability_genetic_score\\t%.4f", cor(liability, genetic_bt)),
		sprintf("cor_BT_genetic_score\\t%.4f", cor(binary, genetic_bt)),
		sprintf("var_QT_explained_by_genetic_score\\t%.4f", cor(quantitative, genetic_qt)^2),
		sprintf("mean_age\\t%.2f", mean(age)),
		sprintf("males\\t%d", sum(sex_code == 1L)),
		sprintf("females\\t%d", sum(sex_code == 2L))
	)
	writeLines(report, "fixture_summary.txt")
	"""

	stub:
	"""
	printf 'FID\\tIID\\tQT\\tBT\\n' > ${params.fixture_prefix}.pheno
	printf 'FID\\tIID\\tSEX\\tAGE\\tPC1\\n' > ${params.fixture_prefix}.covar
	printf 'FID\\tIID\\tAGE\\tPC1\\n' > ${params.fixture_prefix}.qcovar
	printf 'FID\\tIID\\tSEX\\n' > ${params.fixture_prefix}.catcovar
	printf 'samples\\t0\\n' > fixture_summary.txt
	"""
}
