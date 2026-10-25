# Dense LIANA consensus fixture

`liana_dense_consensus.h5ad` exercises dense expression at the production LIANA
rank-aggregate interface. The existing compact CSV panel exceeds LIANA's allowed
missing-resource proportion. This independent synthetic panel covers every
component gene in the pinned consensus resource instead.

## Input contract

- Shape: 40 cells by 1,877 human gene symbols, sorted lexicographically.
- `X`: dense NumPy int32 counts, finite and strictly positive. Baseline counts
  range from 1 to 3. There is no `.raw`, alternative layer or precomputed result.
- `obs.group`: categorical `sender` and `receiver`, with 20 cells each. Use
  `group` as the production process `obs_key`.
- Ordered cells: `cell_000` to `cell_039`. Gene identifiers equal `var.symbols`.
- Five ligands have 24 to 27 counts in every sender. Their receptor components
  have 24 to 27 counts in every receiver. In the opposite group, these genes
  retain baseline counts. All complex subunits are expressed.
- H5AD arrays use gzip level 4 and compatible UTF-8/categorical encodings.

The deliberately enriched consensus interactions are:

| Ligand | Receptor or receptor complex |
| --- | --- |
| CXCL12 | CXCR4 |
| CCL5 | CCR5 |
| TGFB1 | TGFBR1_TGFBR2 |
| IL6 | IL6R_IL6ST |
| MIF | CD74_CXCR4 |

These counts test software behaviour. They are not biological measurements or
evidence of communication. Background expression makes every consensus pair
eligible, including pairs outside this deliberately enriched set.

## Resource provenance

The generator reads `liana/resource/omni_resource.csv` from `liana==1.6.1`, selects
`resource == 'consensus'`, and splits complex names on underscores. The resource
contains 4,624 ligand-receptor pairs and 1,877 unique component genes. All those
genes are expressed, so the missing-resource proportion is zero rather than
near LIANA's 0.98 limit. Only the gene symbols determine the synthetic panel;
the generator does not reuse biological expression data.

The bundled resource matches the [authoritative source at LIANA's 1.6.1 commit](https://github.com/saezlab/liana-py/blob/ea90b7561c549dbd149c0297e28d7ea7d903b55f/src/liana/resource/omni_resource.csv).
Its SHA-256 is
`9b4880f638198725177805d1cf5777b8d2f512444613d579158f1529c196ca2a`.
The generator refuses another LIANA version, a different resource hash or an
unexpected resource size. LIANA's source and bundled resource are distributed
under its [BSD 3-Clause licence](https://github.com/saezlab/liana-py/blob/ea90b7561c549dbd149c0297e28d7ea7d903b55f/LICENSE).

## Reproduction and validation

Run from the repository root into a new directory. The generator also imports
the shared compression validator from `generate_de_raw_fixtures.py`.

```bash
uv run --python 3.12.14 \
    --with liana==1.6.1 --with scanpy==1.12.4 \
    --with anndata==0.12.19 --with numpy==2.4.2 --with pandas==2.3.3 \
    --with scipy==1.17.1 --with h5py==3.16.0 \
    python generate_liana_dense_fixture.py \
    --output-dir /path/to/new-output --validate-liana
```

Generation uses NumPy seed `20261025`. It checks dense positive counts, exact
resource coverage, group sizes, directional enrichment and all complex subunits
before writing. Read-back checks enforce matrix/metadata equality, compression,
string encodings and the 1 MB size limit. Existing output files are never
overwritten. Omitting `--validate-liana` skips only the library inference check.

The optional check runs the real LIANA 1.6.1 `rank_aggregate` on a copy after
`scanpy.pp.log1p`, with `resource_name='consensus'`, `groupby='group'`,
`use_raw=False`, `n_jobs=1`, `n_perms=5` and seed `20261025`. It requires:

- A nonempty `liana_res` with 18,496 rows, representing 4,624 pairs across all
  four source-target group combinations.
- Finite `magnitude_rank` and `specificity_rank` values between zero and one.
- Exactly one sender-to-receiver row for each documented interaction, with
  `lr_means > 3` and `lr_logfc > 1`, matching the deliberately enriched counts.

The checker prints result columns, selected interaction statistics and actual
package versions. The fixture retains its original dense counts and no results.
LIANA can convert its internal working copy to sparse form during inference.
This external generator validation adds no pipeline assertion process. Production
regressions must still execute the real module on the dense input and inspect
its output metadata or tables.

Validation with the pinned command returned a `(18496, 13)` result. Each of the
five directed interactions passed the stated bounds. For example, CXCL12/CXCR4
had `lr_means = 3.266916` and `lr_logfc = 3.225603`. These rounded values describe
the validation run; the robust oracle is the interaction identity and the bounds
above, not exact floating-point equality or a stochastic rank snapshot.
