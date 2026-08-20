# Test data for custom/clustering, clustermetrics and clustervisualization modules

These test files are used by the following custom modules:

- `custom/clustering`
- `custom/clustermetrics`
- `custom/clustervisualization`

The files were generated from the simulated PLINK dataset located in `popgen/plink_simulated.*` (following the same explicit-command convention used in `popgen/README.md`).

**Files and how they were generated:**

- `test.eigenvec`: PCA output (~200 samples) generated with  
  `plink2 --pfile popgen/plink_simulated --pca 10 --out test`

- `test_clusters.csv`: Cluster assignments produced by scikit-learn `KMeans(n_clusters=3, n_init=10, random_state=42)` applied to the PCA feature matrix.

- `test_features.tsv`: Feature matrix (TSV format) consisting of the principal component scores extracted from `test.eigenvec` (without sample identifier columns).  
  **This file is used to calculate the clustering quality metrics** (in the `clustermetrics` module)  
  **and for the cluster visualizations** (in the `clustervisualization` module).

These files test correct input-format parsing and module functionality (technical validity only).
