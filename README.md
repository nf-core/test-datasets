# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the [`nf-core/spatialtranscriptomics`](https://github.com/nf-core/spatialtranscriptomics) pipeline

## Introduction

nf-core is a collection of high quality Nextflow pipelines. This repository contains various files for CI and unit testing of nf-core/spatialtranscriptomics pipeline and infrastructure.

The test data includes a 10x Visium spatial transcriptomics sample [**CID4465**](https://zenodo.org/record/4739739) and a scRNA-seq data from [**GSE176078**: **GSM5354528-CID4465**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078) from the study:

> Wu, S.Z., Al-Eryani, G., Roden, D.L. et al. A single-cell and spatially resolved atlas of human breast cancers. **Nat Genet** 53, 1334–1347 (2021). https://doi.org/10.1038/s41588-021-00911-1

<p align="middle">
	<img src="https://github.com/sdomanskyi/test-datasets/blob/spatialtranscriptomics/testdata/test-dataset/CID4465_ST/spatial/tissue_lowres_image.png?raw=true" width="400"/>
	<img src="https://github.com/sdomanskyi/test-datasets/blob/spatialtranscriptomics/testdata/test-dataset-subsampled/CID4465_tissue_lowres_image_marked.png?raw=true" height="400"/>
	<img src="https://github.com/sdomanskyi/test-datasets/blob/spatialtranscriptomics/testdata/test-dataset-subsampled/CID4465_spots_counts.png?raw=true" height="200"/>
</p>


## Contents

```
testdata
├───test-dataset
│   │   samplesheet.csv
│   │
│   ├───CID4465_SC
│   │   └───filtered_feature_bc_matrix
│   │           barcodes.tsv.gz
│   │           features.tsv.gz
│   │           matrix.mtx.gz
│   │
│   └───CID4465_ST
│       ├───raw_feature_bc_matrix
│       │       barcodes.tsv.gz
│       │       features.tsv.gz
│       │       matrix.mtx.gz
│       │
│       └───spatial
│               aligned_fiducials.jpg
│               detected_tissue_image.jpg
│               scalefactors_json.json
│               tissue_hires_image.png
│               tissue_lowres_image.png
│               tissue_positions_list.csv
│
└───test-dataset-subsampled
    │   ...
    ...
```


## Script written for subsampling

```python
# Script written by Sergii Domanskyi

import os
import gzip
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.io import mmwrite 
import scanpy as sc

def gz(fname):

    '''Compress file with gzip and remove source
    '''

    with open(fname) as f_in:
        with gzip.open(fname + '.gz', 'wt') as f_out:
            f_out.writelines(f_in)

    os.remove(fname)

    return

# Alternatively data could be loaded from mtx of full dataset
dataDir = 'path/to/output/of/full/dataset/'

# Select cell ids from clusters 2 & 4
sc_clusters = pd.read_csv(dataDir + 'STdeconvolve_sc_cluster_ids.csv', index_col=0)
sc_clusters = sc_clusters.loc[(sc_clusters.values==2) | (sc_clusters.values==4)].index

# Load counts data from the full dataset
sc_adata = sc.read(dataDir + 'sc_adata_raw.h5ad')
df_sc = sc_adata.to_df().T[sc_clusters]
df_sc = df_sc.loc[df_sc.index[(df_sc.sum(axis=1)>=80)][:1000]].astype(int)

# Load and mask a subset of spots
st_adata_plain = sc.read(dataDir + 'st_adata_raw.h5ad')
st_adata_sub = st_adata_plain[(st_adata_plain.obsm["spatial"][:, 0] > 6300) 
                            & (st_adata_plain.obsm["spatial"][:, 0] < 7700) 
                            & (st_adata_plain.obsm["spatial"][:, 1] > 7000), :].copy()
sc.pl.spatial(st_adata_sub, img_key="hires", color=["n_counts"])

sc_genes = df_sc.index
sc_and_st_genes = st_adata_sub.var.index.intersection(sc_genes)
st_genes = sc_and_st_genes.union(st_adata_sub.var.loc[~st_adata_sub.var.index.isin(sc_genes)].sample(46, random_state=0).index)
df_st = st_adata_sub.to_df().T.loc[st_genes].astype(int)

saveDataDir = dataDir + 'test-dataset-subsampled/'

# Write scRNA-seq data to mtx format
pd.Series(sc_genes).to_csv(saveDataDir + 'CID4465_SC/filtered_feature_bc_matrix/' + 'features.tsv.gz', sep='\t', index=False, header=False)
pd.Series(df_sc.columns).to_csv(saveDataDir + 'CID4465_SC/filtered_feature_bc_matrix/' + 'barcodes.tsv.gz', sep='\t', index=False, header=False)
fname = saveDataDir + 'CID4465_SC/filtered_feature_bc_matrix/' + 'matrix.mtx'
mmwrite(fname, csr_matrix(df_sc.values), comment='Subsampled')
gz(fname)

# Write ST counts data to mtx format
df_features_sc = st_adata_sub.var.loc[st_genes][['gene_ids', 'feature_types']].reset_index()[['gene_ids', 'index', 'feature_types']]
df_features_sc.to_csv(saveDataDir + 'CID4465_ST/raw_feature_bc_matrix/' + 'features.tsv.gz', sep='\t', index=False, header=False)
pd.Series(df_st.columns).to_csv(saveDataDir + 'CID4465_ST/raw_feature_bc_matrix/' + 'barcodes.tsv.gz', sep='\t', index=False, header=False)
fname = saveDataDir + 'CID4465_ST/raw_feature_bc_matrix/' + 'matrix.mtx'
mmwrite(fname, csr_matrix(df_st.values), comment='metadata_json: {"software_version": "Cell Ranger 4", "format_version": 2}')
gz(fname)

# Exclude from the tissue spots outside selection
df_tisue = pd.read_csv(dataDir + 'test-dataset/CID4465_ST/spatial/tissue_positions_list.csv', index_col=0, header=None)
print(df_tisue.shape[0], df_tisue[1].sum())
df_tisue = df_tisue.loc[df_st.columns]
print(df_tisue.shape[0], df_tisue[1].sum())
df_tisue.reset_index().to_csv(dataDir + 'test-dataset-subsampled/CID4465_ST/spatial/tissue_positions_list.csv', index=False, header=False)

```


## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).

