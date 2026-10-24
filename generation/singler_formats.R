#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(anndataR)
  library(HDF5Array)
  library(Matrix)
  library(SingleCellExperiment)
  library(SingleR)
  library(SummarizedExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 12L) {
  stop("Expected archive, panel, query, output, report, stem, label, stratum, profile and gene options")
}

source_archive <- args[[1L]]
panel_path <- args[[2L]]
query_path <- args[[3L]]
output_dir <- args[[4L]]
report_path <- args[[5L]]
archive_stem <- args[[6L]]
label_field <- args[[7L]]
stratify_field <- args[[8L]]
profile_limit <- as.integer(args[[9L]])
minimum_categories <- as.integer(args[[10L]])
gene_limit <- as.integer(args[[11L]])
panel_hash <- args[[12L]]

stopifnot(file.exists(source_archive), file.exists(panel_path), file.exists(query_path))
stopifnot(profile_limit > 0L, minimum_categories > 0L, gene_limit > 0L)

query_genes <- readLines(panel_path, warn = FALSE)
source_root <- file.path(tempdir(), "singleR-source")
dir.create(source_root, showWarnings = FALSE, recursive = TRUE)
archive_members <- untar(source_archive, list = TRUE)
expected_root <- paste0(archive_stem, "/")
if (!all(startsWith(archive_members, expected_root))) {
  stop("Source tar members do not share the expected archive stem")
}
untar(source_archive, exdir = source_root)
reference_dir <- file.path(source_root, archive_stem)
reference <- loadHDF5SummarizedExperiment(dir = reference_dir)

required_label_fields <- c("label.main", "label.fine")
if (!all(required_label_fields %in% colnames(colData(reference)))) {
  stop("Source reference is missing label.main or label.fine")
}
if (!all(c(label_field, stratify_field) %in% colnames(colData(reference)))) {
  stop("Requested label or profile-stratification field is absent")
}
if (anyNA(colData(reference)[, required_label_fields, drop = FALSE])) {
  stop("Source reference has missing requested label metadata")
}

shared_genes <- query_genes[query_genes %in% rownames(reference)]
if (length(shared_genes) < gene_limit) {
  stop(paste("Only", length(shared_genes), "query genes overlap the reference"))
}
selected_genes <- shared_genes[seq_len(gene_limit)]

strata <- as.character(colData(reference)[[stratify_field]])
profile_count <- min(profile_limit, ncol(reference))
if (profile_count == ncol(reference)) {
  selected_profiles <- seq_len(ncol(reference))
} else {
  categories <- sort(unique(strata))
  by_category <- lapply(categories, function(category) {
    which(strata == category)[order(colnames(reference)[strata == category])]
  })
  names(by_category) <- categories
  selected_profiles <- integer()
  while (length(selected_profiles) < profile_count) {
    progressed <- FALSE
    for (category in categories) {
      candidates <- by_category[[category]]
      if (length(candidates) > 0L && length(selected_profiles) < profile_count) {
        selected_profiles <- c(selected_profiles, candidates[[1L]])
        by_category[[category]] <- candidates[-1L]
        progressed <- TRUE
      }
    }
    if (!progressed) {
      stop("Profile stratification could not select the requested number of profiles")
    }
  }
}

reference_subset <- reference[selected_genes, selected_profiles]
colData(reference_subset) <- colData(reference)[selected_profiles, , drop = FALSE]
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_reference <- file.path(output_dir, archive_stem)
saveHDF5SummarizedExperiment(
  reference_subset,
  dir = output_reference,
  replace = TRUE,
  level = 6L,
  chunkdim = c(min(256L, gene_limit), min(32L, profile_count))
)

reloaded <- loadHDF5SummarizedExperiment(dir = output_reference)
stopifnot(
  identical(as.integer(dim(reloaded)), c(gene_limit, profile_count)),
  all(required_label_fields %in% colnames(colData(reloaded))),
  !anyNA(colData(reloaded)[, required_label_fields, drop = FALSE]),
  length(unique(as.character(colData(reloaded)[[label_field]]))) >= minimum_categories
)

query_adata <- read_h5ad(query_path)
query <- query_adata$as_SingleCellExperiment(x_mapping = "counts", assays_mapping = FALSE)
query_shared_genes <- intersect(rownames(reloaded), rownames(query))
if (length(query_shared_genes) != gene_limit) {
  stop("Reloaded reference genes do not match the query feature names")
}
query <- query[query_shared_genes, , drop = FALSE]
reloaded <- reloaded[query_shared_genes, , drop = FALSE]
assay(reloaded) <- as(assay(reloaded), "dgCMatrix")
predictions <- SingleR(
  test = query,
  assay.type.test = 1L,
  ref = reloaded,
  labels = colData(reloaded)[[label_field]],
  num.threads = 1L
)
if (nrow(predictions) != ncol(query) || all(is.na(predictions$labels))) {
  stop("SingleR did not produce usable query predictions")
}

profile_categories <- length(unique(as.character(colData(reloaded)[[label_field]])))
report <- list(
  archive_stem = archive_stem,
  dimensions = as.integer(dim(reloaded)),
  source_gene_count = length(query_genes),
  shared_gene_count = length(shared_genes),
  gene_count = length(selected_genes),
  ordered_gene_list_sha256 = panel_hash,
  profiles = ncol(reloaded),
  label_field = label_field,
  label_categories = profile_categories,
  stratify_field = stratify_field,
  stratify_categories = length(unique(strata[selected_profiles])),
  query_genes_used = length(query_shared_genes),
  query_cells = ncol(query),
  inference = list(
    query_cells = nrow(predictions),
    predicted_categories = length(unique(as.character(predictions$labels))),
    label_field = label_field
  ),
  runtime = list(
    R = R.version.string,
    SingleR = as.character(packageVersion("SingleR")),
    HDF5Array = as.character(packageVersion("HDF5Array")),
    anndataR = as.character(packageVersion("anndataR"))
  )
)
jsonlite::write_json(report, report_path, auto_unbox = TRUE, pretty = TRUE)
