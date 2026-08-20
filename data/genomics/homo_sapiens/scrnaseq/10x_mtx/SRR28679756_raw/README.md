## Gene expression matrix in 10X Matrix Market format

Conversion of modules/data/genomics/homo_sapiens/scrnaseq/h5ad/SRR28679756_raw_matrix.h5ad into 10X-style Matrix Market format, i.e. 3 files:

- matrix.mtx.gz: Matrix Market sparse matrix
- barcodes.tsv.gz: 1 column, no header, with cell barcodes that are the column names of matrix.mtx.gz
- features.tsv.gz: 3 columns, no header, with gene_name, gene_name, "Gene Expression"

Conventionally the first column of features.tsv.gz contains gene IDs, but SRR28679756_raw_matrix.h5ad only contains gene names.
