"""Generate and validate matched expression-source regression fixtures."""

from __future__ import annotations

import argparse
import platform
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, issparse

SEED = 20261025
N_CELLS = 64
N_GENES = 128
N_MARKERS = 24
FILENAMES = ("de_raw_absent.h5ad", "de_raw_conflicting.h5ad")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def make_pair() -> tuple[ad.AnnData, ad.AnnData]:
    rng = np.random.default_rng(SEED)
    counts = rng.integers(0, 4, size=(N_CELLS, N_GENES), dtype=np.int32)
    counts[:, : 2 * N_MARKERS] = 0
    for offset, rows in ((0, slice(0, 32)), (N_MARKERS, slice(32, 64))):
        counts[rows, offset : offset + N_MARKERS] = rng.integers(
            8, 12, size=(32, N_MARKERS), dtype=np.int32
        )
        counts[rows, offset] += 12

    genes = (
        [f"A_MARKER_{i:02d}" for i in range(N_MARKERS)]
        + [f"B_MARKER_{i:02d}" for i in range(N_MARKERS)]
        + [f"BACKGROUND_{i:03d}" for i in range(N_GENES - 2 * N_MARKERS)]
    )
    obs = pd.DataFrame(
        {"group": pd.Categorical(["A"] * 32 + ["B"] * 32, categories=["A", "B"])},
        index=pd.Index([f"cell_{i:03d}" for i in range(N_CELLS)], name="cell_id"),
    )
    var = pd.DataFrame(
        {"symbols": np.asarray(genes, dtype=object)},
        index=pd.Index(genes, name="gene_id"),
    )
    absent = ad.AnnData(X=csr_matrix(counts), obs=obs, var=var)
    conflicting = absent.copy()
    reversed_counts = counts.copy()
    reversed_counts[:, :N_MARKERS] = counts[:, N_MARKERS : 2 * N_MARKERS]
    reversed_counts[:, N_MARKERS : 2 * N_MARKERS] = counts[:, :N_MARKERS]
    # Raw holds log-expression, as expected by Scanpy when it selects raw implicitly.
    conflicting.raw = ad.AnnData(
        X=csr_matrix(np.log1p(reversed_counts.astype(np.float32))),
        obs=obs.copy(),
        var=var.copy(),
    )
    return absent, conflicting


def validate_pair(absent: ad.AnnData, conflicting: ad.AnnData) -> None:
    require(absent.shape == conflicting.shape == (N_CELLS, N_GENES), "Unexpected shape")
    require(
        absent.raw is None and conflicting.raw is not None, "Unexpected raw presence"
    )
    require(absent.obs_names.is_unique and absent.var_names.is_unique, "Duplicate axes")
    pd.testing.assert_frame_equal(absent.obs, conflicting.obs)
    pd.testing.assert_frame_equal(absent.var, conflicting.var)
    require(
        absent.obs["group"].value_counts().to_dict() == {"A": 32, "B": 32},
        "Groups differ",
    )
    require(issparse(absent.X) and issparse(conflicting.X), "X must remain sparse")
    require((absent.X != conflicting.X).nnz == 0, "X differs between fixtures")
    require(absent.X.dtype == conflicting.X.dtype == np.int32, "Counts must be int32")
    require(
        np.isfinite(absent.X.data).all() and (absent.X.data >= 0).all(),
        "Invalid counts",
    )
    for data in (absent, conflicting):
        require(
            not data.layers and not data.uns,
            "Unexpected layers or unstructured metadata",
        )
        require(
            not data.obsm and not data.varm and not data.obsp and not data.varp,
            "Unexpected slots",
        )

    raw = conflicting.raw
    require(issparse(raw.X) and raw.X.dtype == np.float32, "Raw must be sparse float32")
    pd.testing.assert_index_equal(raw.obs_names, absent.obs_names)
    pd.testing.assert_frame_equal(raw.var, absent.var)
    x = absent.X.toarray()
    expected = x.copy()
    expected[:, :N_MARKERS] = x[:, N_MARKERS : 2 * N_MARKERS]
    expected[:, N_MARKERS : 2 * N_MARKERS] = x[:, :N_MARKERS]
    np.testing.assert_array_equal(
        raw.X.toarray(), np.log1p(expected.astype(np.float32))
    )

    for group, rows, other_rows, offset in (
        ("A", slice(0, 32), slice(32, 64), 0),
        ("B", slice(32, 64), slice(0, 32), N_MARKERS),
    ):
        columns = slice(offset, offset + N_MARKERS)
        require((x[rows, columns] >= 8).all(), f"{group} markers lack signal in X")
        require((x[other_rows, columns] == 0).all(), f"{group} markers leak in X")
        require(
            (raw.X[rows, columns].toarray() == 0).all(),
            f"{group} raw signal not reversed",
        )
        require(
            (raw.X[other_rows, columns].toarray() > 0).all(),
            f"{group} raw markers absent",
        )
        require(
            (
                x[rows, offset] > x[rows, offset + 1 : offset + N_MARKERS].max(axis=1)
            ).all(),
            f"{group} primary marker is not strongest",
        )


def validate_hdf5(path: Path) -> None:
    with h5py.File(path, "r") as handle:

        def visit(name: str, obj: h5py.Dataset | h5py.Group) -> None:
            if isinstance(obj, h5py.Dataset) and obj.shape != () and obj.size:
                require(
                    obj.compression == "gzip",
                    f"{path.name}:{name} is not gzip compressed",
                )
            encoding = obj.attrs.get("encoding-type", "")
            require(
                encoding != "nullable-string-array",
                f"{path.name}:{name} uses nullable strings",
            )
            if encoding == "string-array":
                require(
                    h5py.check_string_dtype(obj.dtype).encoding == "utf-8",
                    f"{name} is not UTF-8",
                )

        handle.visititems(visit)
    require(path.stat().st_size < 1_000_000, f"{path.name} exceeds the 1 MB target")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    paths = [args.output_dir / name for name in FILENAMES]
    for path in paths:
        if path.exists():
            raise FileExistsError(f"Refusing to overwrite {path}")
    pair = make_pair()
    validate_pair(*pair)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for data, path in zip(pair, paths, strict=True):
        data.write_h5ad(path, compression="gzip", compression_opts=4)
        validate_hdf5(path)
    validate_pair(*(ad.read_h5ad(path) for path in paths))
    print(f"Validated matched {N_CELLS}-cell/{N_GENES}-gene pair; seed={SEED}")
    print(f"python=={platform.python_version()}")
    for package in ("anndata", "numpy", "pandas", "scipy", "h5py"):
        print(f"{package}=={version(package)}")
    for path in paths:
        print(f"{path}: {path.stat().st_size} bytes")


if __name__ == "__main__":
    main()
