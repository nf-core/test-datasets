#!/usr/bin/env Rscript

library(Seurat)

sce <- readRDS("pbmc3k_processed_sce.rds")

seurat_obj <- as.Seurat(sce, data = NULL)

saveRDS(seurat_obj, file = "pbmc3k_processed_seurat.rds")
