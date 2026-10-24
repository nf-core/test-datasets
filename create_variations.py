#!/usr/bin/env python3
"""Build small, independent AnnData fixtures for structural tests."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any, cast

os.environ.setdefault("ANNDATA_ALLOW_WRITE_NULLABLE_STRINGS", "0")

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix, csr_matrix

N_CELLS = 128
N_GENES = 128
LABELS = np.array(["T-Cell", "B cell", "CD4+ T", "CD8+ T", "NK cell", "Monocyte"])
CONDITIONS = np.array(["Healthy", "Disease", "Treated"])


def _obs(adata: ad.AnnData) -> pd.DataFrame:
    return cast(pd.DataFrame, adata.obs)


def _var(adata: ad.AnnData) -> pd.DataFrame:
    return cast(pd.DataFrame, adata.var)


def _counts(seed: int, case_number: int) -> tuple[np.ndarray, np.random.Generator]:
    rng = np.random.default_rng(np.random.SeedSequence([seed, case_number]))
    values = rng.poisson(2.0, size=(N_CELLS, N_GENES)).astype(np.int32)
    values[rng.random(values.shape) < 0.82] = 0
    values[np.arange(N_CELLS), np.arange(N_CELLS)] += 1
    return values, rng


def _new_adata(
    seed: int, case_number: int
) -> tuple[ad.AnnData, np.ndarray, np.random.Generator]:
    counts, rng = _counts(seed, case_number)
    obs = pd.DataFrame(
        {
            "sample": np.array(
                [f"sample_{i % 2}" for i in range(N_CELLS)], dtype=object
            ),
        },
        index=pd.Index([f"cell_{i:03d}" for i in range(N_CELLS)]),
    )
    genes = [f"GENE{i:03d}" for i in range(N_GENES)]
    var = pd.DataFrame(
        {
            "symbols": np.array(genes, dtype=object),
            "gene_id": np.array(
                [f"ENSG{i:011d}" for i in range(N_GENES)], dtype=object
            ),
        },
        index=pd.Index(genes),
    )
    return ad.AnnData(X=csr_matrix(counts), obs=obs, var=var), counts, rng


def _with_layers(adata: ad.AnnData, counts: np.ndarray) -> None:
    adata.X = csr_matrix(counts.astype(np.float32) * 0.5 + 0.25)
    adata.layers["counts"] = csr_matrix(counts.copy())
    adata.layers["spliced"] = csr_matrix(counts + 2)
    adata.layers["unspliced"] = csr_matrix(counts * 2 + 3)


def _add_optional_slots(adata: ad.AnnData, rng: np.random.Generator) -> None:
    adata.obsm["X_pca"] = rng.normal(size=(N_CELLS, 4)).astype(np.float32)
    adata.obsm["X_umap"] = rng.normal(size=(N_CELLS, 2)).astype(np.float32)
    adata.varm["PCs"] = rng.normal(size=(N_GENES, 4)).astype(np.float32)
    adata.uns["neighbors"] = {"connectivities_key": "connectivities"}
    adata.uns["pca"] = {"variance_ratio": rng.random(4).astype(np.float32)}
    rows = np.arange(N_CELLS - 1)
    cols = rows + 1
    values = np.ones(N_CELLS - 1, dtype=np.float32)
    adata.obsp["connectivities"] = csr_matrix(
        (
            np.concatenate([values, values]),
            (np.concatenate([rows, cols]), np.concatenate([cols, rows])),
        ),
        shape=(N_CELLS, N_CELLS),
    )


def _fixture_cases(seed: int) -> tuple[dict[str, ad.AnnData], dict[str, str]]:
    cases: dict[str, ad.AnnData] = {}
    distinctions: dict[str, str] = {}

    def add(name: str, description: str, make: Any) -> None:
        case_number = len(cases) + 1
        adata, counts, rng = _new_adata(seed, case_number)
        make(adata, counts, rng)
        cases[f"{name}.h5ad"] = adata
        distinctions[f"{name}.h5ad"] = description

    add(
        "counts_in_X",
        "Counts are stored in X; there are no layers.",
        lambda a, _x, _r: None,
    )

    def counts_in_layer(
        adata: ad.AnnData, counts: np.ndarray, _rng: np.random.Generator
    ) -> None:
        adata.X = csr_matrix(counts.astype(np.float32) * 0.25 + 0.125)
        adata.layers["counts"] = csr_matrix(counts.copy())

    add(
        "counts_in_layer",
        "X differs from the selected integer counts layer.",
        counts_in_layer,
    )

    def symbols_as_index(
        adata: ad.AnnData, _counts: np.ndarray, _rng: np.random.Generator
    ) -> None:
        _var(adata).index = _var(adata)["symbols"].astype(str)
        adata.var = _var(adata).drop(columns="symbols")

    add(
        "symbols_as_index",
        "Gene symbols are represented only by var_names.",
        symbols_as_index,
    )
    add(
        "symbols_in_column",
        "Gene symbols are in the non-boolean var column symbols.",
        lambda _a, _x, _r: None,
    )

    def duplicates(
        adata: ad.AnnData, counts: np.ndarray, _rng: np.random.Generator
    ) -> None:
        symbols = _var(adata)["symbols"].to_numpy(dtype=object)
        for i in range(8):
            symbols[64 + i] = symbols[i]
            counts[i, i] = i + 1
            counts[i, 64 + i] = i + 10
        _var(adata)["symbols"] = symbols
        adata.X = csr_matrix(counts)

    add(
        "with_duplicates",
        "Eight pairs collide after selecting var.symbols; their count columns differ.",
        duplicates,
    )

    def isoforms(
        adata: ad.AnnData, counts: np.ndarray, _rng: np.random.Generator
    ) -> None:
        names = _var(adata).index.to_numpy(dtype=object)
        for i in range(8):
            names[i] = f"GENE{i:03d}.1"
            names[64 + i] = f"GENE{i:03d}.2"
            counts[i, i] = i + 1
            counts[i, 64 + i] = i + 10
        _var(adata).index = names
        adata.X = csr_matrix(counts)

    add(
        "with_isoforms",
        "Eight numeric-suffix pairs collide after isoform suffix removal.",
        isoforms,
    )

    def batch_different(
        adata: ad.AnnData, _counts: np.ndarray, _rng: np.random.Generator
    ) -> None:
        _obs(adata)["batch_id"] = _obs(adata)["sample"].copy()
        adata.obs = _obs(adata).drop(columns="sample")

    add(
        "batch_different_name",
        "Batch labels are present only in obs.batch_id.",
        batch_different,
    )

    def batch_correct(
        adata: ad.AnnData, _counts: np.ndarray, _rng: np.random.Generator
    ) -> None:
        _obs(adata)["batch"] = _obs(adata)["sample"].copy()
        adata.obs = _obs(adata).drop(columns="sample")

    add(
        "batch_correct_name",
        "Batch labels are present only in obs.batch.",
        batch_correct,
    )
    add(
        "batch_missing",
        "Neither sample nor batch labels are present.",
        lambda a, _x, _r: setattr(a, "obs", _obs(a).drop(columns="sample")),
    )

    def batch_conflict(
        adata: ad.AnnData, _counts: np.ndarray, _rng: np.random.Generator
    ) -> None:
        adata.obs["batch"] = "conflict_batch"

    add(
        "batch_conflict",
        "Both sample and batch are present with conflicting values.",
        batch_conflict,
    )

    def label_different(
        adata: ad.AnnData, _counts: np.ndarray, rng: np.random.Generator
    ) -> None:
        _obs(adata)["cell_type"] = rng.choice(LABELS, size=N_CELLS).astype(object)

    add(
        "label_different_name",
        "Labels are present only in obs.cell_type.",
        label_different,
    )

    def label_correct(
        adata: ad.AnnData, _counts: np.ndarray, rng: np.random.Generator
    ) -> None:
        _obs(adata)["label"] = rng.choice(LABELS, size=N_CELLS).astype(object)

    add("label_correct_name", "Labels are present only in obs.label.", label_correct)
    add("label_missing", "No label column is present.", lambda _a, _x, _r: None)

    def label_unknown(
        adata: ad.AnnData, _counts: np.ndarray, rng: np.random.Generator
    ) -> None:
        labels = rng.choice(LABELS, size=N_CELLS).astype(object)
        labels[:16] = "Unassigned"
        _obs(adata)["cell_type"] = labels

    add(
        "label_with_unknown",
        "cell_type includes an explicit Unassigned value.",
        label_unknown,
    )

    def label_nan(
        adata: ad.AnnData, _counts: np.ndarray, rng: np.random.Generator
    ) -> None:
        labels = rng.choice(LABELS, size=N_CELLS).astype(object)
        labels[:12] = np.nan
        _obs(adata)["cell_type"] = labels

    add(
        "label_with_nan",
        "cell_type contains both labels and missing values.",
        label_nan,
    )

    def label_special(
        adata: ad.AnnData, _counts: np.ndarray, _rng: np.random.Generator
    ) -> None:
        _obs(adata)["cell_type"] = np.resize(LABELS, N_CELLS).astype(object)

    add(
        "label_special_chars",
        "cell_type includes spaces, plus signs and hyphens.",
        label_special,
    )

    def label_conflict(
        adata: ad.AnnData, _counts: np.ndarray, rng: np.random.Generator
    ) -> None:
        _obs(adata)["label"] = "conflict_label"
        _obs(adata)["cell_type"] = rng.choice(LABELS, size=N_CELLS).astype(object)

    add("label_conflict", "Both label and cell_type are present.", label_conflict)

    def condition_different(
        adata: ad.AnnData, _counts: np.ndarray, rng: np.random.Generator
    ) -> None:
        _obs(adata)["disease_state"] = rng.choice(CONDITIONS, size=N_CELLS).astype(
            object
        )

    add(
        "condition_different_name",
        "Conditions are present only in obs.disease_state.",
        condition_different,
    )

    def condition_correct(
        adata: ad.AnnData, _counts: np.ndarray, rng: np.random.Generator
    ) -> None:
        _obs(adata)["condition"] = rng.choice(CONDITIONS, size=N_CELLS).astype(object)

    add(
        "condition_correct_name",
        "Conditions are present only in obs.condition.",
        condition_correct,
    )
    add("condition_missing", "No condition column is present.", lambda _a, _x, _r: None)

    def condition_conflict(
        adata: ad.AnnData, _counts: np.ndarray, rng: np.random.Generator
    ) -> None:
        _obs(adata)["condition"] = "conflict_condition"
        _obs(adata)["disease_state"] = rng.choice(CONDITIONS, size=N_CELLS).astype(
            object
        )

    add(
        "condition_conflict",
        "Both condition and disease_state are present.",
        condition_conflict,
    )

    def minimal(
        adata: ad.AnnData, _counts: np.ndarray, _rng: np.random.Generator
    ) -> None:
        adata.obs = _obs(adata).drop(columns="sample")
        adata.layers.clear()

    add(
        "minimal_structure",
        "Only X, obs and var are present; no sample or optional slots.",
        minimal,
    )

    def matrix_csr(
        adata: ad.AnnData, counts: np.ndarray, _rng: np.random.Generator
    ) -> None:
        adata.X = csr_matrix(counts)

    add("matrix_csr", "X is encoded as CSR sparse data.", matrix_csr)

    def matrix_csc(
        adata: ad.AnnData, counts: np.ndarray, _rng: np.random.Generator
    ) -> None:
        adata.X = csc_matrix(counts)

    add("matrix_csc", "X is encoded as CSC sparse data.", matrix_csc)

    def duplicate_cells(
        adata: ad.AnnData, _counts: np.ndarray, _rng: np.random.Generator
    ) -> None:
        names = _obs(adata).index.to_list()
        for i in range(8):
            names[64 + i] = names[i]
        adata.obs_names = names

    add("obs_duplicate_names", "Eight cell-name pairs are duplicated.", duplicate_cells)

    def sample_different(
        adata: ad.AnnData, _counts: np.ndarray, _rng: np.random.Generator
    ) -> None:
        _obs(adata)["sample"] = "different_sample_id"

    add(
        "sample_exists_different",
        "sample is present and differs from the sample identifier.",
        sample_different,
    )
    add(
        "sample_missing",
        "The sample column is absent.",
        lambda a, _x, _r: setattr(a, "obs", _obs(a).drop(columns="sample")),
    )

    def with_obsm(
        adata: ad.AnnData, counts: np.ndarray, rng: np.random.Generator
    ) -> None:
        _with_layers(adata, counts)
        _add_optional_slots(adata, rng)

    add(
        "with_obsm",
        "obsm, varm, uns, obsp and unrelated layers are populated for removal checks.",
        with_obsm,
    )

    def with_layers(
        adata: ad.AnnData, counts: np.ndarray, _rng: np.random.Generator
    ) -> None:
        adata.layers.clear()
        adata.X = csr_matrix(counts.astype(np.float32) * 0.5 + 0.25)
        adata.layers["counts"] = csr_matrix(counts.copy())
        adata.layers["spliced"] = csr_matrix(counts + 2)
        adata.layers["unspliced"] = csr_matrix(counts * 2 + 3)

    add(
        "with_layers",
        "X, counts, spliced and unspliced contain intentionally distinct matrices.",
        with_layers,
    )

    def complex_structure(
        adata: ad.AnnData, counts: np.ndarray, rng: np.random.Generator
    ) -> None:
        _with_layers(adata, counts)
        _add_optional_slots(adata, rng)
        _obs(adata)["batch_id"] = _obs(adata)["sample"].copy()
        _obs(adata)["cell_type"] = rng.choice(LABELS, size=N_CELLS).astype(object)
        _var(adata)["feature_type"] = np.array(
            ["Gene Expression"] * N_GENES, dtype=object
        )

    add(
        "complex_structure",
        "Counts layers, string var metadata and every optional AnnData slot are populated.",
        complex_structure,
    )

    if len({id(adata) for adata in cases.values()}) != len(cases):
        raise AssertionError(
            "Each structural fixture must use an independent AnnData object"
        )
    return cases, distinctions


def _check_compression(path: Path) -> dict[str, int]:
    checked = 0
    exempt = 0
    with h5py.File(path, "r") as handle:

        def visit(_name: str, obj: h5py.Dataset | h5py.Group) -> None:
            nonlocal checked, exempt
            if not isinstance(obj, h5py.Dataset):
                return
            if obj.shape == () or obj.size == 0:
                exempt += 1
                return
            checked += 1
            if obj.compression != "gzip":
                raise AssertionError(
                    f"{path.name}: dataset {_name} uses {obj.compression!r}, not gzip"
                )

        handle.visititems(visit)
    return {"gzip_datasets": checked, "scalar_or_empty_exemptions": exempt}


def _matrix_equal(left: Any, right: Any) -> bool:
    difference = csr_matrix(left) - csr_matrix(right)
    return difference.nnz == 0


def _columns_equal(matrix: Any, left: int, right: int) -> bool:
    return _matrix_equal(matrix[:, left], matrix[:, right])


def build_variations(
    output_dir: str | Path,
    seed: int = 42,
    *,
    writer: Any = None,
    write_manifest: bool = True,
) -> dict[str, Any]:
    """Write structural fixtures and return their deterministic QA manifest."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    cases, distinctions = _fixture_cases(seed)
    manifest: dict[str, Any] = {
        "seed": seed,
        "shape": [N_CELLS, N_GENES],
        "synthetic_counts": True,
        "raw_counts_layer_written": False,
        "files": {},
    }

    if writer is None:

        def writer(adata: ad.AnnData, path: Path) -> None:
            adata.write_h5ad(path, compression="gzip", compression_opts=4)

    for filename, adata in cases.items():
        path = output_dir / filename
        manifest_path = output_dir / "structural_manifest.json"
        if path.exists() or manifest_path.exists():
            raise FileExistsError(
                f"Refusing to overwrite existing output: {path if path.exists() else manifest_path}"
            )
        writer(adata, path)
        if path.stat().st_size >= 5_000_000:
            raise AssertionError(f"{filename} is not below the 5 MB hard limit")
        if path.stat().st_size >= 1_000_000:
            raise AssertionError(f"{filename} exceeds the 1 MB target")

        restored = ad.read_h5ad(path)
        compression = _check_compression(path)
        file_qa: dict[str, Any] = {
            "bytes": path.stat().st_size,
            "shape": list(restored.shape),
            "x_encoding": type(restored.X).__name__,
            "obs_columns": list(restored.obs.columns),
            "var_columns": list(restored.var.columns),
            "layers": list(restored.layers.keys()),
            "optional_slots": {
                "obsm": list(restored.obsm.keys()),
                "varm": list(restored.varm.keys()),
                "obsp": list(restored.obsp.keys()),
                "uns": list(restored.uns.keys()),
            },
            "obs_names_unique": bool(restored.obs_names.is_unique),
            "var_names_unique": bool(restored.var_names.is_unique),
            "symbols_dtype": (
                str(restored.var["symbols"].dtype)
                if "symbols" in restored.var
                else None
            ),
            "compression": compression,
            "distinction": distinctions[filename],
        }
        if restored.shape != (N_CELLS, N_GENES):
            raise AssertionError(f"{filename} has unexpected shape {restored.shape}")
        manifest["files"][filename] = file_qa

    files = manifest["files"]
    if files["matrix_csr.h5ad"]["x_encoding"] != "csr_matrix":
        raise AssertionError("matrix_csr.h5ad did not round-trip as CSR")
    if files["matrix_csc.h5ad"]["x_encoding"] != "csc_matrix":
        raise AssertionError("matrix_csc.h5ad did not round-trip as CSC")
    if files["with_duplicates.h5ad"]["var_names_unique"] is False:
        raise AssertionError(
            "Duplicate-symbol fixture must have unique source var_names"
        )
    if "symbols" not in files["symbols_in_column.h5ad"]["var_columns"]:
        raise AssertionError("symbols_in_column.h5ad is missing var.symbols")
    if files["symbols_in_column.h5ad"]["symbols_dtype"] in {"bool", "boolean"}:
        raise AssertionError(
            "var.symbols must remain non-boolean for dtype-validation tests"
        )
    if files["obs_duplicate_names.h5ad"]["obs_names_unique"]:
        raise AssertionError(
            "obs_duplicate_names.h5ad does not contain duplicate cell names"
        )
    if files["with_obsm.h5ad"]["optional_slots"]["obsm"] == []:
        raise AssertionError("with_obsm.h5ad must populate optional slots")

    counts_layer = ad.read_h5ad(output_dir / "counts_in_layer.h5ad")
    if _matrix_equal(counts_layer.X, counts_layer.layers["counts"]):
        raise AssertionError("counts_in_layer.h5ad must distinguish X and counts")
    layered = ad.read_h5ad(output_dir / "with_layers.h5ad")
    if not all(
        not _matrix_equal(left, right)
        for left, right in (
            (layered.X, layered.layers["counts"]),
            (layered.layers["counts"], layered.layers["spliced"]),
            (layered.layers["spliced"], layered.layers["unspliced"]),
        )
    ):
        raise AssertionError(
            "with_layers.h5ad must contain distinct X/counts/spliced/unspliced matrices"
        )

    for filename, adata in cases.items():
        if filename == "with_duplicates.h5ad":
            duplicate_count = int(_var(adata)["symbols"].duplicated().sum())
            if duplicate_count != 8:
                raise AssertionError(
                    f"Expected eight duplicate symbols, found {duplicate_count}"
                )
            distinct_pairs = sum(
                not _columns_equal(adata.X, i, 64 + i) for i in range(8)
            )
            if distinct_pairs != 8:
                raise AssertionError(
                    f"Expected eight distinct duplicate-gene count pairs, found {distinct_pairs}"
                )
            files[filename]["duplicate_symbols"] = duplicate_count
            files[filename]["distinct_duplicate_gene_count_pairs"] = distinct_pairs
        if filename == "with_isoforms.h5ad":
            collapsed = pd.Series(adata.var_names).str.replace(r"\.\d+", "", regex=True)
            collisions = int(collapsed.duplicated().sum())
            if collisions != 8:
                raise AssertionError(
                    f"Expected eight isoform collisions, found {collisions}"
                )
            distinct_pairs = sum(
                not _columns_equal(adata.X, i, 64 + i) for i in range(8)
            )
            if distinct_pairs != 8:
                raise AssertionError(
                    f"Expected eight distinct isoform count pairs, found {distinct_pairs}"
                )
            files[filename]["isoform_collisions_after_suffix_removal"] = collisions
            files[filename]["distinct_isoform_count_pairs"] = distinct_pairs
        if filename == "label_with_nan.h5ad":
            missing = int(_obs(adata)["cell_type"].isna().sum())
            if missing != 12:
                raise AssertionError(f"Expected twelve missing labels, found {missing}")
            files[filename]["missing_labels"] = missing

    manifest_path = output_dir / "structural_manifest.json"
    if write_manifest:
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--seed", default=42, type=int)
    args = parser.parse_args()
    manifest = build_variations(args.output_dir, seed=args.seed)
    largest = max(
        (details["bytes"], name) for name, details in manifest["files"].items()
    )
    print(f"Built {len(manifest['files'])} structural fixtures in {args.output_dir}")
    print(f"Largest fixture: {largest[1]} ({largest[0]} bytes)")
    print(f"Manifest: {args.output_dir / 'structural_manifest.json'}")


if __name__ == "__main__":
    main()
