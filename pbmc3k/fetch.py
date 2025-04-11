#!/usr/bin/env python

import scanpy as sc

adata_raw = sc.datasets.pbmc3k()
adata_raw.write_h5ad('pbmc3k_raw.h5ad')
adata_raw.var.index.name = None

adata_raw.var["symbols"] = adata_raw.var_names
adata_raw.var.index = adata_raw.var["gene_ids"]
adata_raw.var.index.name = None
adata_raw.write_h5ad('pbmc3k_raw_geneid.h5ad')

adata_processed = sc.datasets.pbmc3k_processed()
adata_processed.var.index.name = None
adata_processed.write_h5ad('pbmc3k_processed.h5ad')
sc.pp.sample(adata_processed, n=500, copy=True, rng=0).to_df().to_csv('pbmc3k_processed.csv')

geneid_mapping = sc.datasets.pbmc3k().var["gene_ids"].to_dict()
adata_processed.var.index = adata_processed.var.index.map(geneid_mapping)
adata_processed.write_h5ad('pbmc3k_processed_geneid.h5ad')
