# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)
Test data to be used for automated testing with the [`nf-core/spatialtranscriptomics`](https://github.com/nf-core/spatialtranscriptomics) pipeline

## Introduction

nf-core is a collection of high quality Nextflow pipelines. This repository contains various files for CI and unit testing of nf-core/spatialtranscriptomics pipeline and infrastructure.

The test data includes a 10x Visium spatial transcriptomics sample [**CID4465**](https://zenodo.org/record/4739739) and a scRNA-seq data from [**GSE176078**: **GSM5354528-CID4465**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078) from the study:

> Wu, S.Z., Al-Eryani, G., Roden, D.L. et al. A single-cell and spatially resolved atlas of human breast cancers. **Nat Genet** 53, 1334–1347 (2021). https://doi.org/10.1038/s41588-021-00911-1

The subsampled CID4465_ST dataset contains 187 spots with 1000 genes, of which 126 are in tissue. The subsampled CID4465_SC dataset contains 273 cells with 1000 genes. The intersection of CID4465_ST genes and CID4465_SC genes is a set of 954 genes. For details on subsampling see [script](https://github.com/sdomanskyi/test-datasets/blob/spatialtranscriptomics/testdata/test-dataset-subsampled/subsample.py).

<p align="middle">
	<img src="https://github.com/sdomanskyi/test-datasets/blob/spatialtranscriptomics/testdata/test-dataset/CID4465_ST/spatial/tissue_lowres_image.png?raw=true" width="300"/>
	<img src="https://github.com/sdomanskyi/test-datasets/blob/spatialtranscriptomics/testdata/test-dataset-subsampled/CID4465_tissue_lowres_image_marked.png?raw=true" height="300"/>
	<img src="https://github.com/sdomanskyi/test-datasets/blob/spatialtranscriptomics/testdata/test-dataset-subsampled/CID4465_spots_counts.png?raw=true" height="175"/>
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


## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).

