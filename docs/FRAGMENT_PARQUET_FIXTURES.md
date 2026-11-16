# Fragment Parquet fixtures

These supplemental fixtures support the `feat/parquet-fragments` pipeline migration.
They do not replace the legacy `.pkl` files or change the main generation stages.
The package contains 36 Parquet tables, two CSV inputs and a machine-readable
[manifest](../intermediates/fragments/parquet_manifest.json).
Every fixture is smaller than 1,000,000 bytes.

## Generate and check

Run the supplemental generator after the existing dataset generation stages:

```bash
uv run generate_fragment_parquet_fixtures.py --dataset-root .
uv run generate_fragment_parquet_fixtures.py --dataset-root . --check
```

The script declares Python 3.12, NumPy 2.3.5, pandas 3.0.0 and PyArrow 23.0.1 in
its inline dependency metadata. `uv` manages the environment without installing
system packages. pandas 3 is needed to load the existing categorical Leiden
pickle; pandas 2.3 cannot restore its categorical state. The repository's
`requirements.txt` already includes PyArrow and remains unchanged.

Generation reads the trusted local legacy pickles and the first column of
`pbmc/SRR28679756_filtered_matrix.csv`. It does not download data, run models,
modify source fixtures or invoke the pipeline. Re-running it replaces only its
own outputs. `--dataset-root` selects an existing dataset checkout, not an empty
output directory.

`--check` writes nothing. It rebuilds the expected tables in memory, checks their
values and physical schema, then compares every stored file and the manifest
byte for byte. Both modes serialise each table twice to check deterministic
bytes in the pinned environment. Validation checks pandas metadata, column
order, identifier validity, categorical values, category order, unused
categories, nullable types, null masks and Snappy compression. The manifest
records each relative path, row count, ordered fields, Arrow and pandas types,
null counts, categorical metadata, byte size, SHA-256, source and consumer tests.
It also records input SHA-256 values and the pipeline revision from which the
synthetic payloads were copied: `a219db8f144eb1f8be17104b8b4b733539396298`.

## Physical representation

All Parquet tables store `__index__` as the first ordinary column. They are
written with `frame.rename_axis('__index__').reset_index().to_parquet(...)`,
using `index=False`, `engine='pyarrow'` and `compression='snappy'`. There is no
hidden pandas index. The generator uses the equivalent assignment to
`frame.index.name` before `reset_index()`.

`obs`, `var`, `obsm`, prediction and count-table identifiers are unique, non-null
strings, except the explicitly malformed reader fixtures. Arrow stores these
strings as `large_string` in the pinned environment. Embedding columns have
string names `0` through `N-1`, with their original `float32` values. Categories
retain pandas metadata, their complete category lists and the ordered flag.
pandas 3 may restore an object-backed category index as a string-backed index;
category values, codes, unused levels and ordering remain unchanged.

Tabular `uns` retains its original numeric index in `__index__`. The synthetic
LIANA-like table uses the deliberately non-contiguous index `[7, 3]`; it must not
be replaced with cell IDs or a fresh range. Nested objects, including PAGA, RGG
and decoupler results, stay pickled. The current legacy fixture inventory
contains the nested `scvi-global-1.0_paga.pkl`, which is retained without a
Parquet counterpart. All existing pickles remain available to older pipelines.

## Legacy counterparts

Paths in this table are relative to `intermediates/`.

| Path | Rows | Columns after `__index__` |
| --- | ---: | --- |
| `fragments/X_scvi.parquet` | 16 | `0` to `7`: float32 |
| `fragments/X_scvi-global_umap.parquet` | 12 | `0`, `1`: float32 |
| `fragments/scvi-global-1.0_entropy.parquet` | 12 | `scvi-global-1.0:entropy`: float32 |
| `fragments/scvi-global-1.0_leiden.parquet` | 12 | `scvi-global-1.0_leiden`: unordered category `['0', '1']` |
| `doublet_detection/test_doubletdetection.parquet` | 384 | `test_doubletdetection`: bool |
| `doublet_detection/test_scrublet.parquet` | 384 | `test_scrublet`: bool |
| `doublet_detection/test_solo.parquet` | 384 | `test_solo`: bool |

Each counterpart preserves source values, row order and column order. Only the
explicit index representation and embedding column-name types change.

## Regression fixtures

The following paths are relative to `intermediates/fragments/regressions/`.
Strings use Arrow `large_string`, float64 uses `double`, float32 uses `float`,
integers use `int64` and booleans use `bool`. `__index__` is a string column unless
a row explicitly describes an intentional error.

| File | Rows | Columns after `__index__` |
| --- | ---: | --- |
| `a_obs.parquet` | 3 | `label`: string, `n_counts`, `score`: float64, `flag`: nullable bool, `version.1`: string |
| `b_obs.parquet` | 2 | `label`: string, `n_counts`: float64, `score`: int64, `flag`: bool, `version.1`: string |
| `a_var.parquet` | 1 | `selected`: bool, `score`: float64 |
| `b_var.parquet` | 2 | `selected`: bool, `score`: float64 |
| `a_conflict.parquet`, `b_conflict.parquet` | 1 each | `label`: string, respectively `T_cell` and `B_cell` |
| `a_score.parquet`, `b_score.parquet` | 1 each | `score`: int64, respectively `1` and `2` |
| `duplicate_index.parquet` | 2 | `label`: string; both identifiers equal the first base cell |
| `missing_index.parquet` | 1 | `label`: string; `__index__` deliberately absent |
| `null_index.parquet` | 2 | `label`: string; second identifier deliberately null |
| `annotations.parquet` | 2 | `annotation:test:per_cell`: string, both `T_cell` |
| `missing_batch.parquet` | 1 | `regression_batch`: string, `batch_a` |
| `empty_batch.parquet` | 1 | `empty_batch`: float64, null, for `not_a_cell` |
| `labels.parquet` | 128 | `safe`, `collision`, `path`, `traversal`, `staged_input`: strings |
| `panels.parquet` | 128 | `panel_1`, `panel_10`, `panel_19`, `panel_20`, `panel_32`: bool |
| `panels_plot.parquet` | 32 | Same panel columns, restricted to the first 32 genes |
| `previous_qc.parquet` | 128 | `pct_counts_in_top_20_genes`: float64, all `17.0` |
| `positive_counts.parquet` | 8 | `GENE000` to `GENE031`: int64 |
| `zero_predictions.parquet` | 128 | `first`, `second`: bool, all false |
| `small_predictions.parquet` | 128 | `first`, `second`: bool, true only for `cell_000` |
| `unpaired_obs.parquet` | 12 | `unpaired_donor`: string |
| `interleaved.parquet`, `contiguous.parquet` | 64 each | `regression_batch`: string |
| `groups.parquet` | 64 | `interaction_group`: string |
| `categorical_obs.parquet`, `categorical_var.parquet` | 3 each | Category and nullable columns described below |
| `large_integer.parquet` | 2 | `large_integer`: nullable pandas `Int64`, values `[9007199254740993, null]` |

`large_integer.parquet` uses the first two base cell IDs in their original
order. Its integer value is `2**53 + 1`, which cannot round-trip through
float64 without losing precision. The physical `large_integer` field is Arrow
`int64`, with one null and pandas nullable `Int64` metadata.

`labels.parquet` retains all dangerous labels, punctuation and case exactly as
the original split-column test. The safe labels cycle through `A B`, `a_b`,
`T:CD4+` and `B-cell.(naive)`. The other columns retain normalisation collisions,
`bad/name`, `..` and the staged-input name `labels`.

The panel masks select the first 1, 10, 19, 20 or 32 genes in original order.
The positive matrix uses `1 + (cell_number + gene_number) % 7` and cell names
`cell_0` through `cell_7`. The accompanying `positive_counts.csv` contains the
same values, a `cell` header and the original trailing newline. Use the CSV
for `ADATA_READCSV`; the Parquet counterpart documents the exact table without
requiring a new expression-input reader.

`duplicate_columns.csv` preserves the original invalid `cell,label,label`
header and its single `T_cell,B_cell` row. pandas cannot write duplicate column
names to Parquet. Use this CSV only for legacy CSV-header rejection; use the
missing-index and null-index Parquet fixtures for the new reader rejection
tests.

The donor mapping preserves loop order: donor `D01`, `D02`, `D03`, then label
`TypeA`, `TypeB`, then condition `Healthy`, `Disease`. IDs have the form
`D01__TypeA__Healthy`, with value `D01_Healthy`.

Scanorama and LIANA fixtures read cell IDs from the dense source CSV without
sorting. Interleaved labels alternate `z_batch`, `a_batch`; contiguous labels
use 32 `z_batch` cells followed by 32 `a_batch` cells. Interaction groups
alternate `sender`, `receiver` in the same row order.

Both categorical fixtures contain:

- `ordered_label`: `['high', null, 'low']`, ordered categories
  `['low', 'medium', 'high', 'unused']`.
- `unordered_label`: `['B', 'A', null]`, unordered categories
  `['B', 'A', 'unused']`.
- `nullable_integer`: `[1, null, 3]`, pandas `Int64`.
- `nullable_boolean`: `[true, null, false]`, pandas `boolean`.
- `nullable_string`: `['001', null, 'NA']`, pandas `string`.
- `nullable_float`: `[1.5, null, 0.0]`, float64.

The observation IDs are the first two base cells followed by `not_a_cell`.
The variable IDs are `FAM177A1P1`, `POMP`, `not_a_gene`.

`intermediates/fragments/regressions/uns_table.parquet` has two rows and columns
`__index__`: int64, `source`, `target`, `ligand_complex`, `receptor_complex`:
strings, `magnitude_rank`, `specificity_rank`: float64. Its values are
`(7, sender, receiver, CXCL12, CXCR4, 0.125, 0.25)` and
`(3, receiver, sender, CCL5, CCR5, 0.5, null)`. These are synthetic structural
values, not inferred interactions.

## Doublet score and merge regressions

`generate_doublet_parquet_regressions.py` adds 28 tables in
`intermediates/doublet_detection/regressions/` for the regression tests introduced
in scdownstream PRs 393 and 395. Run it with:

```bash
uv run generate_doublet_parquet_regressions.py
```

The generator checks the table round trips and Snappy compression. Its
`manifest.json` records paths, row counts, columns and SHA-256 values.
These fixtures use 128 `cell_000` to `cell_127` identifiers. Calls retain the
original reversed or interleaved row order. Score tables retain reordered
native score/class columns and a missing Scrublet score. `sample_a` and
`sample_b` tables use the sample-prefixed identifiers from `ADATA_UNIFY`.
Malformed tables cover duplicate, missing and extra identifiers, multiple
prediction columns, null calls, ambiguous strings and numeric calls.
The existing fixtures and generation stages remain unchanged.

## Pipeline consumers

The manifest maps fixtures to these 13 test files. It records the source test
for each synthetic payload; generation does not depend on a pipeline checkout.

| Test path | Fixtures |
| --- | --- |
| `modules/local/adata/extend/tests/main.nf.test` | Compatible, conflicting, malformed and categorical fragments; legacy counterparts; tabular `uns` |
| `modules/local/adata/entropy/tests/main.nf.test` | `empty_batch` |
| `modules/local/adata/splitcol/tests/labels.nf.test` | `labels` |
| `modules/local/custom/aggregate_per_cell_annotations/tests/main.nf.test` | `annotations` |
| `modules/local/custom/doubletremoval/tests/main.nf.test` | Three legacy prediction counterparts |
| `modules/local/custom/doubletremoval/tests/plot_edge_cases.nf.test` | `zero_predictions`, `small_predictions` |
| `modules/local/edgepython/differential/tests/design.nf.test` | `unpaired_obs` |
| `modules/local/liana/rankaggregate/tests/dense.nf.test` | `groups` |
| `modules/local/scanorama/integrate/tests/main.nf.test` | `interleaved`, `contiguous` |
| `modules/local/scanorama/integrate/tests/batch_validation.nf.test` | `missing_batch` |
| `modules/local/scanpy/filter/tests/small_panels.nf.test` | `panels`, `previous_qc` |
| `modules/local/scanpy/filter/tests/small_panels_plot.nf.test` | `positive_counts`, `panels_plot` |
| `subworkflows/local/finalize/tests/main.nf.test` | Legacy counterparts and tabular `uns` |

## Local pipeline testing before publication

The published fork base URL is
`https://raw.githubusercontent.com/nictru/test-datasets/scdownstream`.
The root samplesheets and `singleR/references.csv` use this same fork for
embedded dataset URLs. The SingleR generator copies `singleR/references.csv`,
so its generated sheet retains the fork URLs. Historical commit-pinned source
URLs remain unchanged.

Use the dataset checkout as the test-data base, without a trailing slash:

```text
--pipelines_testdata_base_path /nfs/proj/n.trummer/opencode/worktrees/test-datasets-fragment-parquet
```

Nextflow `file()` accepts local absolute paths, so no HTTP server or fork push
is required. Run tests from the pipeline checkout and retain their usual
container profile. The shared NFS checkout must be visible to execution nodes.

If a consumer specifically requires HTTP, serve this checkout read-only:

```bash
python3 -m http.server 8765 --bind 127.0.0.1 \
    --directory /nfs/proj/n.trummer/opencode/worktrees/test-datasets-fragment-parquet
```

Then use `http://127.0.0.1:8765` as the base URL for local execution. Loopback
does not reach another worker node. Prefer the absolute NFS path for cluster
tests. The generator never commits files or publishes repository changes.
