# Matched differential-expression fixtures

`de_raw_absent.h5ad` and `de_raw_conflicting.h5ad` test implicit selection of
AnnData `.raw` in Scanpy ranking, marker filtering and dot plots. Both files
contain the same `X`, ordered axes and metadata. Only the second contains `.raw`.
These synthetic counts do not represent biological measurements.

## Contract

- Shape: 64 cells by 128 genes. `obs.group` contains ordered categories `A` and
  `B`, with 32 cells each. Use `group` as the ranking `obs_key`.
- Cells are `cell_000` to `cell_063`. Gene identifiers equal `var.symbols`.
- `X` contains non-negative int32 CSR counts. Each group expresses 24 exclusive
  markers at 8 to 11 counts. `A_MARKER_00` and `B_MARKER_00` receive 12 additional
  counts and are unambiguous primary markers of their respective groups.
- The other 80 genes contain seeded background counts from 0 to 3. Every cell
  has a positive total. The panel supports real ranking and the default marker
  plots, including PCA before the dot plot.
- Conflicting `.raw.X` contains float32 natural-log `log1p` expression after
  exchanging the A and B marker blocks. Its axes and gene metadata match `X`.
  Background expression is unchanged apart from `log1p`. The primary marker
  identities therefore reverse when Scanpy selects `.raw` implicitly.
- Neither file contains layers, embeddings, graphs or precomputed rankings.
  H5AD arrays use gzip level 4. Strings use compatible UTF-8/string-array and
  categorical encodings, not nullable-string-array.

The production ranking module applies `log1p` to `X`. Explicit `use_raw=False`
must give matching rankings, filtered markers and plots for these files.
The 24 exclusive markers per group pass its in-group/out-group fraction filters.
This contract applies to Wilcoxon and supported t-tests; it does not supply
missing statistics for logistic regression.
Perfectly separated markers can have tied Wilcoxon scores. Compare marker sets
or gene-keyed table rows rather than assuming a unique first-ranked marker.

## Reproduction

Run from the repository root into a new directory. No global installation is
required. The validated environment uses CPython 3.12.14 and the package versions
in this command:

```bash
uv run --python 3.12.14 \
    --with anndata==0.12.19 --with numpy==2.4.2 --with pandas==2.3.3 \
    --with scipy==1.17.1 --with h5py==3.16.0 \
    python generate_de_raw_fixtures.py --output-dir /path/to/new-output
```

The fixed NumPy seed is `20261025`. The generator validates the pair before
writing and after reading it back. It checks exact count equality, metadata and
ordered axes, group sizes, marker separation, the complete reversed raw matrix,
gzip compression, string encodings and the 1 MB per-file size limit. Existing
outputs are never overwritten. Validation reads expression arrays only within
this external fixture recipe; pipeline regressions inspect output metadata and
differential-expression tables.

Reproduction guarantees matrix and metadata equality with the pinned environment.
HDF5 container bytes can differ between library builds. The generator prints
the actual runtime versions and file sizes.

An additional local compatibility check used Scanpy 1.12.3 and AnnData 0.12.10.
Wilcoxon and t-test produced identical explicit-X ranking tables and 24 filtered
markers per group for the matched inputs. Rank plots and dot plots succeeded.
Implicit raw selection reversed the A marker set in both methods. Scanpy emitted
expected warnings about tied-score plot limits and omitted two-group dendrograms.
This check exercised library calls, not the complete pipeline process.

Publish both files and this generator together. Pipeline regressions must use
`https://raw.githubusercontent.com/nictru/test-datasets/<commit>/anndata-variations/<filename>`
with an immutable fixture commit. After upstream merges the data, switch those
two references to the immutable upstream commit. Other test-data URLs retain
their existing base.
