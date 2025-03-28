# Test Data

Produced data was created from subsetting the files found in the test data of the [ARACNe3](https://github.com/califano-lab/ARACNe3/tree/3d8791a23e3bd8fd0d74f3b8d48f912e81d00f14/data) github repo, to have smaller data and shorter test run times.

```console
califano-lab/ARACNe3
├── ...
└── data
    ├── exp_mat.txt.gz
    ├── regulators.txt
    └── standardized_results_1.txt
```

## Regulator Gene List

List of regulator genes for which the regulon network will be created. `/data/regulators.txt` was generated from taking the first 100 genes from the tool's repo [`regulators.txt`](https://github.com/califano-lab/ARACNe3/blob/3d8791a23e3bd8fd0d74f3b8d48f912e81d00f14/data/regulators.txt) file.

## Expression Data Matrix

Gene (rows) expression TPMs per sample (column). `/data/expression_matrix.tsv` was generated from taking 5 first samples (columns) from the tool's repo [`exp_mat.txt.gz`](https://github.com/califano-lab/ARACNe3/blob/3d8791a23e3bd8fd0d74f3b8d48f912e81d00f14/data/exp_mat.txt.gz) file.

