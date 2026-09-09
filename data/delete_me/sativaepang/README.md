# sativaepang test fixtures

Two tiers, for `sativaepang/reference`, `sativaepang/lootasks`, `sativaepang/looplace` and
`sativaepang/looscore` module tests (nf-core/modules), and the `sativa` subworkflow that
chains them with `raxmlng/search`.

## Tiny: sativa-epang's own bundled example

`sativaepang_tiny.phy` (PHYLIP), `sativaepang_tiny.tax`, `sativaepang_tiny_synonym.txt`:
the 38-sequence example bundled with
[`Aaramis/sativa-epang`](https://github.com/Aaramis/sativa-epang/tree/main/example)
(`example/test.phy`/`test.tax`/`synonym.txt`, unmodified), also used by the tool's own
bioconda recipe test. Fast smoke test; the tool's own default config reports 0 mislabels
on this set.

`sativaepang_tiny.raxml.bestTree`/`sativaepang_tiny.raxml.bestModel`: RAxML-NG 2.0.3
(`raxml-ng --search --msa sativaepang_tiny.phy --model GTR+G --seed 42`) run against the
alignment above, for `sativaepang/reference`'s `-reftree`/`-refmodel` inputs.

`sativaepang_tiny.refjson`/`sativaepang_tiny.model`: `sativaepang/reference`'s own output
(`sativa-epang -s sativaepang_tiny.phy -t sativaepang_tiny.tax -x bac -reftree
sativaepang_tiny.raxml.bestTree -refmodel sativaepang_tiny.raxml.bestModel -stage
reference`), for `sativaepang/lootasks`'s input.

## Small: GTDB archaeal 16S dataset, with a known injected mislabel pair

`gtdb_archaea_16s_aligned.fasta`/`gtdb_archaea_16s.tax`: the same 121-sequence,
10-species GTDB archaeal 16S dataset documented on the `sativa` branch's own README
(`## GTDB archaeal 16S dataset`), copied unmodified. Crucially, it already carries a
deliberate positive control: `DupHaloA`/`DupSulfoA`, a swapped-label pair maximally
divergent in the tree (`Haloferax volcanii` content tagged as `Saccharolobus
islandicus`'s lineage and vice versa). Confirmed 2026-09-08: the full `sativa-epang`
leave-one-out chain (below) flags both at Phylum level with high confidence
(0.851/0.819), plus three of the genuine hard cases the `sativa` branch README already
documents as honest raxtax disagreements (`Haloarcula`/`Haloferax` family confusion,
one `Haloferax` species pair) — a real, meaningful mislabel-detection assertion for
tests, not just "the process completed."

`gtdb_archaea_16s.raxml.bestTree`/`gtdb_archaea_16s.raxml.bestModel`: RAxML-NG 2.0.3
(`raxml-ng --search --msa gtdb_archaea_16s_aligned.fasta --model GTR+G --seed 42`), same
role as the tiny tier's. RAxML-NG warns of 71 near-zero branches on this alignment
(short, gappy 16S fragments over a 10-species panel) — expected given the dataset's own
purpose (a taxonomic diversity ladder, not a densely-sampled clade), not a defect in the
tree used for placement testing.

`gtdb_archaea_16s.refjson`/`gtdb_archaea_16s.model`: `sativaepang/reference`'s output on
this dataset, same invocation pattern as the tiny tier's (`-x bac`), for
`sativaepang/lootasks`'s input.

## A known sativa-epang bug this fixture set exposed

`sativa-epang -stage loo-tasks -r <refjson>` (no `SATIVA_EPANG_MODEL` set) embeds the
`.model` file's raw content — including RAxML-NG's own trailing `, name = range`
partition clause — verbatim into `manifest.json`, which `-stage loo-place` then feeds to
EPA-ng's inline `-m` mode per fold. Inline mode rejects that trailing clause
(`Wrong model specification`), while file mode (used at `-stage reference`) requires it
— there is no single format that satisfies both ends of the round trip. All four module
tests work around this by exporting `SATIVA_EPANG_MODEL` themselves, trimmed to the bare
model expression before the first comma (`cut -d ',' -f1 model_file`). Reported to
Auguste (repo owner) directly; not yet fixed upstream as of 2026-09-08.
