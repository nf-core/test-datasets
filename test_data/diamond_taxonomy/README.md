# Test DIAMOND taxonomy databases for nf-core/metatdenovo

Two tiny, self-contained DIAMOND databases used to exercise the pipeline's `--diamond_dbs`
option (see the [nf-core/metatdenovo usage docs](https://nf-co.re/metatdenovo/usage#taxonomic-annotation-with-diamond))
in `nf-test`, without needing the full-size NCBI RefSeq or GTDB databases (which are far too
large to download/store/run in CI).

- `ncbi-refseq-test.dmnd` + `ncbi_taxdump/` — mimics `--diamond_dbs` usage with a real NCBI
  Taxonomy-numbered database and `parse_with_taxdump: true`.
- `gtdb-r220-test.dmnd` + `gtdb_taxdump/` — mimics usage with a GTDB-style database and an
  explicit `ranks` list.

## How these were built (2026-07-25)

nf-core/metatdenovo's `-profile test` (megahit + prodigal) was run once against the two real,
full-size databases (NCBI RefSeq protein, snapshot 2024-01-15, and GTDB R09-RS220) on a separate
machine, giving a ground-truth mapping of every predicted ORF (`megahit.prodigal.faa.gz`) to a
real taxid for each database (`diamond blastp --outfmt 102` + `taxonkit lineage`).

From that ground truth, a handful of taxa were picked per database — the most abundant hits plus
a few rarer ones, spanning both prokaryotes and (for NCBI) a eukaryote — and, for each, up to 5
of the **actual ORF protein sequences** (≥30 aa) that produced that ground-truth hit were pulled
out of the same `megahit.prodigal.faa.gz` and used directly as the reference sequences in the
tiny test database, tagged with the real taxid from the full-database run. No sequences were
downloaded from NCBI or GTDB — every reference sequence here is exactly one of the ORFs that
`nf-core/metatdenovo`'s own `tests/default.nf.test` (megahit + prodigal, `-profile test`)
produces, so the provenance of every sequence is fully reproducible by just running that test.

Taxa included:

| Database  | Taxon                            | Rank        |
| --------- | --------------------------------- | ----------- |
| NCBI      | *Saccharomyces cerevisiae* S288C  | species     |
| NCBI      | *Saccharomyces paradoxus*         | species     |
| NCBI      | *Bacteroides fragilis*            | species     |
| NCBI      | *Acinetobacter baumannii*         | species     |
| NCBI      | *Escherichia coli*                | species     |
| NCBI      | *Pseudomonas aeruginosa*          | species     |
| NCBI      | *Klebsiella pneumoniae*           | species     |
| GTDB      | *Acinetobacter*                   | genus       |
| GTDB      | Enterobacteriaceae                | family      |
| GTDB      | *Bacteroides*                     | genus       |
| GTDB      | *Escherichia*                     | genus       |
| GTDB      | *Acinetobacter pittii* (000369045)| species     |
| GTDB      | *Parabacteroides intestinipullorum* (019114965) | species |
| GTDB      | *Citrobacter*                     | genus       |

`ncbi_taxdump/{names,nodes}.dmp` is a verbatim subset (just the lines needed for these taxa and
their full ancestor lineage up to root) of a real NCBI taxdump — real NCBI taxids throughout,
so e.g. taxid `562` still means *E. coli* as it does everywhere else.

`gtdb_taxdump/{names,nodes}.dmp` uses **surrogate small-integer taxids** (root=1, others 2-22)
instead of GTDB's real numeric ids (which run up into the billions, e.g. `2091866253`). This
isn't for size (see caveat below) — GTDB's real numeric ids are already synthetic/arbitrary in
GTDB itself (the taxon *name* is what's meaningful there, not the number), so nothing scientific
is lost; the taxon names in `names.dmp` are still the real ones.

## Known limitation

Because each reference sequence is a distinct real ORF mapped 1:1 to a single taxid (not an
engineered near-duplicate across taxa), no query should score comparably against references from
two different taxa. This means DIAMOND's own internal tie-triggered LCA collapse (the mechanism
that, against the real full-size databases, produces higher-rank assignments like
"Enterobacteriaceae" when several species tie for best hit) is **not** exercised by this fixture.
That's a gap in what this specific test can verify, but it's DIAMOND's own algorithm rather than
pipeline code — what this fixture does correctly exercise is the pipeline's own logic: the
`diamond blastp` / `taxonkit lineage` invocation, summary-table generation and `ranks` parsing,
`parse_with_taxdump`, and handling of lineage strings of varying depth (several of the chosen
GTDB taxa are themselves genus/family-rank, not species, which does test that).

## Regenerating or extending

Run `nf-core/metatdenovo`'s `tests/default.nf.test` to reproduce `megahit.prodigal.faa.gz`, then
pick ORF sequences for whichever taxa are needed and tag them with taxids (real, for NCBI-style
databases; a small surrogate range, for GTDB-style ones) via `diamond makedb --taxonmap`, using a
subset of a real taxdump for `--taxonnames`/`--taxonnodes` covering the chosen taxa's full
lineage up to root.
