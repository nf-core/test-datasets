*.contrasts.yaml : These are files created from the original csv contrasts files converted into YAML format.
MaxQuant_contrasts.yaml created from https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/proteomics/maxquant/MaxQuant_contrasts.csv
SRP254919.contrasts.yaml created from https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/mus_musculus/rnaseq_expression/SRP254919.contrasts.csv

SRP254919.propd_grea_test_sets.gmt : A small Ensembl-id GMT with two synthetic gene sets drawn from the first 50 gene IDs of SRP254919.salmon.merged.gene_counts.top1000cov.tsv (in the `modules` branch). Used by the nf-core/differentialabundance `test_rnaseq_propd_grea` profile to exercise propr/grea functional enrichment, which operates directly on the abundance matrix's Ensembl IDs and so cannot reuse the symbol-based mh.all.v2022.1.Mm.symbols.gmt.
