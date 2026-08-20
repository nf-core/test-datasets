# orf_catalogue

Minimal crafted ORF catalogue used as input to the `custom/orfcollapse` module
test. The real chr20 data does not contain the case orfcollapse handles (the
same micropeptide encoded at distinct, non-overlapping loci), so this small
fixture is hand-built to exercise it: two identical-peptide small ORFs on
opposite strands/loci, one canonical CDS and one unique small ORF.

`orf_class` is positional only — it records where an ORF sits relative to the
annotated CDS and never encodes length. The `is_smorf` column carries the
length flag (`aa_length` <= `--smorf-max-aa`), and `orf_type_native` preserves
each caller's own ORF-type label. The two identical-peptide rows share a class
deliberately: orfcollapse prefers the more specific class when a peptide cluster
spans several, so differing classes would change which row survives the fold.

- `cohort.catalogue.{bed12,tsv}`, `cohort.orf_to_gene.tsv`, `cohort.catalogue.mqc.tsv`: orfmerge-shaped catalogue.
- `cohort.catalogue.aa.fasta`: catalogue peptides (headers carry the `(+)/(-)` suffix bedtools getfasta -nameOnly -s emits).
- `cohort_cluster.tsv`: MMseqs2 cluster TSV grouping the two identical small ORFs.
