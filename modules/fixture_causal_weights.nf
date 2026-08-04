// Draw the causal variants and their effect sizes for the simulated traits.
//
// Two disjoint causal sets are drawn from the merged fixture variants, one for the
// quantitative trait and one for the binary trait's genetic liability, so the two
// traits are not collinear. Every draw goes through a fixed seed with the RNG kind
// pinned explicitly, so the weight file is byte-identical on every run.
process FIXTURE_CAUSAL_WEIGHTS {
	container "community.wave.seqera.io/library/r-base:4.4.3--1e564c44feffeaa0"
	publishDir params.outdir_fixture_pheno_cov, mode: 'copy'

	input:
	path pvar

	output:
	path "causal_weights.tsv", emit: weights

	script:
	"""
	#!/usr/bin/env Rscript
	set.seed(${params.fixture_seed},
			 kind = "Mersenne-Twister",
			 normal.kind = "Inversion",
			 sample.kind = "Rejection")

	variants <- read.table("${pvar}", comment.char = "#", header = FALSE,
						   stringsAsFactors = FALSE)
	colnames(variants) <- c("CHROM", "POS", "ID", "REF", "ALT")
	n_variants <- nrow(variants)

	n_causal <- ${params.fixture_n_causal}
	if (2 * n_causal > n_variants) {
		stop("not enough variants for two disjoint causal sets")
	}

	drawn <- sample.int(n_variants, 2 * n_causal)
	causal_qt <- sort(drawn[seq_len(n_causal)])
	causal_bt <- sort(drawn[n_causal + seq_len(n_causal)])

	weights <- data.frame(ID = variants[["ID"]],
						  A1 = variants[["ALT"]],
						  BETA_QT = 0,
						  BETA_BT = 0,
						  stringsAsFactors = FALSE)
	weights[["BETA_QT"]][causal_qt] <- round(rnorm(n_causal), 6)
	weights[["BETA_BT"]][causal_bt] <- round(rnorm(n_causal), 6)

	write.table(weights, file = "causal_weights.tsv", sep = "\\t", quote = FALSE,
				row.names = FALSE, col.names = TRUE)
	"""

	stub:
	"""
	printf 'ID\\tA1\\tBETA_QT\\tBETA_BT\\n' > causal_weights.tsv
	"""
}
