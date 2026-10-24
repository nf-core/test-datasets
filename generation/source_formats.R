#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 10) {
    stop("Usage: source_formats.R FORMAT FILTERED RAW FILTERED_OUT RAW_OUT SAMPLE NCELLS NGENES NBACKGROUND GENE_PANEL")
}

format <- args[[1]]
filtered_input <- args[[2]]
raw_input <- args[[3]]
filtered_output <- args[[4]]
raw_output <- args[[5]]
sample_id <- args[[6]]
n_cells <- as.integer(args[[7]])
n_genes <- as.integer(args[[8]])
n_background <- as.integer(args[[9]])
gene_panel <- args[[10]]

suppressPackageStartupMessages({
    library(Matrix)
    library(SingleCellExperiment)
    library(Seurat)
})

counts_for <- function(object, format) {
    if (format == "sce") {
        assay_name <- if ("counts" %in% SummarizedExperiment::assayNames(object)) "counts" else SummarizedExperiment::assayNames(object)[[1]]
        return(SummarizedExperiment::assay(object, assay_name))
    }
    Seurat::GetAssayData(object, assay = Seurat::DefaultAssay(object), layer = "counts")
}

subset_object <- function(object, format, genes, cells) {
    if (format == "sce") {
        return(object[genes, cells])
    }
    object[genes, cells]
}

selected_cells <- function(counts, n) {
    depth <- Matrix::colSums(counts)
    order(-depth, colnames(counts), method = "radix")[seq_len(n)]
}

selected_genes <- function(counts, n) {
    required <- c("ANXA1", "JCHAIN", "GPNMB", "HSPA12B", "VCAN")
    means <- Matrix::rowMeans(counts)
    means_sq <- Matrix::rowMeans(counts ^ 2)
    ranked <- order(-(means_sq - means ^ 2), rownames(counts), method = "radix")
    required_indices <- match(required, rownames(counts), nomatch = 0L)
    chosen <- unique(c(required_indices[required_indices > 0], ranked))
    chosen[seq_len(n)]
}

design <- function(barcodes) {
    groups <- expand.grid(
        donor = sprintf("D%02d", 1:3),
        condition = c("Healthy", "Disease"),
        label = c("TypeA", "TypeB"),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    groups[((seq_along(barcodes) - 1L) %% nrow(groups)) + 1L, , drop = FALSE]
}

annotate <- function(object, format, sample_id, selected, selected_design = NULL) {
    barcodes <- colnames(object)
    data <- if (format == "sce") as.data.frame(SummarizedExperiment::colData(object)) else object[[]]
    canonical <- c("batch", "donor", "condition", "label")
    data <- data[, setdiff(colnames(data), canonical), drop = FALSE]
    data$source_sample <- sample_id
    data$source_barcode <- barcodes
    data$sample <- sample_id
    data$is_selected_filtered_cell <- selected
    data$source_role <- ifelse(selected, "filtered_cell", "raw_background")
    assigned <- design(barcodes[selected])
    for (column in c("donor", "condition", "label")) {
        alias <- paste0("synthetic_", column)
        data[[alias]] <- "background"
        data[selected, alias] <- if (is.null(selected_design)) assigned[[column]] else selected_design[barcodes[selected], alias]
    }
    counts <- counts_for(object, format)
    data$original_total_counts <- as.integer(Matrix::colSums(counts))
    data$original_n_genes_by_counts <- as.integer(Matrix::colSums(counts > 0))
    if (format == "sce") {
        SummarizedExperiment::rowData(object)$symbols <- rownames(object)
        SummarizedExperiment::colData(object) <- S4Vectors::DataFrame(data, row.names = barcodes)
    } else {
        for (column in canonical) {
            object[[column]] <- NULL
        }
        assay_name <- Seurat::DefaultAssay(object)
        object[[assay_name]] <- Seurat::AddMetaData(
            object[[assay_name]], metadata = rownames(object), col.name = "symbols"
        )
        object <- Seurat::AddMetaData(object, metadata = data)
    }
    object
}

filtered <- readRDS(filtered_input)
raw <- readRDS(raw_input)
filtered_counts <- counts_for(filtered, format)
raw_counts <- counts_for(raw, format)
if (ncol(filtered_counts) < n_cells) {
    stop("Filtered source has too few cells")
}
filtered_indices <- selected_cells(filtered_counts, n_cells)
genes <- readLines(gene_panel, warn = FALSE)
if (length(genes) != n_genes || anyDuplicated(genes) || !all(genes %in% rownames(filtered_counts))) {
    stop("Selected H5AD gene panel is not an authentic compatible RDS feature axis")
}
selected_barcodes <- colnames(filtered_counts)[filtered_indices]
raw_selected_indices <- match(selected_barcodes, colnames(raw_counts), nomatch = 0L)
if (any(raw_selected_indices == 0L)) {
    stop("Raw source lacks selected filtered barcodes")
}
raw_depth <- Matrix::colSums(raw_counts[genes, , drop = FALSE])
background <- which(!(colnames(raw_counts) %in% colnames(filtered_counts)) & raw_depth > 0 & raw_depth < 100)
if (length(background) < n_background) {
    stop("Raw source lacks enough non-zero background droplets below 100 counts")
}
shallow <- background[raw_depth[background] == 2L]
if (length(shallow) < 256L) {
    stop("Raw source lacks 256 two-count droplets on the selected gene panel")
}
shallow <- shallow[order(raw_depth[shallow], colnames(raw_counts)[shallow], method = "radix")][seq_len(256L)]
remaining <- background[!(background %in% shallow)]
deep <- remaining[order(-raw_depth[remaining], colnames(raw_counts)[remaining], method = "radix")][seq_len(n_background - 256L)]
background <- c(shallow, deep)

filtered_out <- subset_object(filtered, format, genes, selected_barcodes)
filtered_out <- annotate(filtered_out, format, sample_id, rep(TRUE, n_cells))
filtered_design <- if (format == "sce") as.data.frame(SummarizedExperiment::colData(filtered_out)) else filtered_out[[]]
rownames(filtered_design) <- filtered_design$source_barcode
raw_barcodes <- colnames(raw_counts)[c(raw_selected_indices, background)]
raw_out <- subset_object(raw, format, genes, raw_barcodes)
raw_out <- annotate(
    raw_out,
    format,
    sample_id,
    seq_along(raw_barcodes) <= length(raw_selected_indices),
    filtered_design[, c("synthetic_donor", "synthetic_condition", "synthetic_label"), drop = FALSE]
)

metadata_for <- function(object, format) {
    if (format == "sce") as.data.frame(SummarizedExperiment::colData(object)) else object[[]]
}

for (object in list(filtered_out, raw_out)) {
    present <- intersect(c("batch", "donor", "condition", "label"), colnames(metadata_for(object, format)))
    if (length(present)) stop(sprintf("Individual source retains canonical metadata: %s", paste(present, collapse = ", ")))
}
raw_metadata <- metadata_for(raw_out, format)
raw_selected <- raw_metadata[raw_metadata$is_selected_filtered_cell, , drop = FALSE]
matched <- filtered_design[match(raw_selected$source_barcode, rownames(filtered_design)), , drop = FALSE]
for (column in c("synthetic_donor", "synthetic_condition", "synthetic_label")) {
    if (!identical(as.character(raw_selected[[column]]), as.character(matched[[column]]))) {
        stop(sprintf("Raw selected cells lost %s", column))
    }
}

saveRDS(filtered_out, filtered_output, compress = "xz")
saveRDS(raw_out, raw_output, compress = "xz")
