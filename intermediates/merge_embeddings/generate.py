#!/usr/bin/env python3
"""Generate deterministic, tiny fixtures for ADATA_MERGEEMBEDDINGS."""

from __future__ import annotations

import argparse
from collections.abc import Callable
from pathlib import Path
from typing import Any, cast

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

N_BASE_CELLS = 50
N_INTEGRATED_CELLS = 30
N_GENES = 100
EMBEDDING_DIM = 10
FILENAMES = [
    f"{kind}_{algorithm}.h5ad"
    for algorithm in ("scvi", "scanvi")
    for kind in ("base", "integrated", "combined")
]


def _obs(adata: ad.AnnData) -> pd.DataFrame:
    return cast(pd.DataFrame, adata.obs)


def create_base_adata(
    rng: np.random.Generator,
    integration_name: str = "scvi",
    include_scanvi_labels: bool = False,
) -> ad.AnnData:
    """Create base data with X_<integration> and unique observation names."""
    cell_names = [f"base_cell_{i}" for i in range(N_BASE_CELLS)]
    gene_names = [f"gene_{i}" for i in range(N_GENES)]
    counts = csr_matrix(rng.poisson(5, size=(N_BASE_CELLS, N_GENES)).astype(np.int32))
    adata = ad.AnnData(
        X=counts,
        obs=pd.DataFrame({"base_only": "base"}, index=pd.Index(cell_names)),
        var=pd.DataFrame(index=pd.Index(gene_names)),
    )
    adata.obsm[f"X_{integration_name}"] = rng.normal(
        size=(N_BASE_CELLS, EMBEDDING_DIM)
    ).astype(np.float32)
    if include_scanvi_labels:
        labels = rng.choice(["TypeA", "TypeB", "TypeC"], size=N_BASE_CELLS)
        adata.obs["label:scANVI"] = pd.Categorical(labels)
    return adata


def create_integrated_adata(
    rng: np.random.Generator, include_scanvi_labels: bool = False
) -> ad.AnnData:
    """Create integrated data with X_emb and unique observation names."""
    cell_names = [f"integrated_cell_{i}" for i in range(N_INTEGRATED_CELLS)]
    gene_names = [f"gene_{i}" for i in range(N_GENES)]
    counts = csr_matrix(
        rng.poisson(5, size=(N_INTEGRATED_CELLS, N_GENES)).astype(np.int32)
    )
    adata = ad.AnnData(
        X=counts,
        obs=pd.DataFrame({"integrated_only": "integrated"}, index=pd.Index(cell_names)),
        var=pd.DataFrame(index=pd.Index(gene_names)),
    )
    adata.obsm["X_emb"] = rng.normal(size=(N_INTEGRATED_CELLS, EMBEDDING_DIM)).astype(
        np.float32
    )
    if include_scanvi_labels:
        labels = rng.choice(["TypeA", "TypeB", "TypeC"], size=N_INTEGRATED_CELLS)
        adata.obs["label:scANVI"] = pd.Categorical(labels)
    return adata


def create_combined_adata(
    base_adata: ad.AnnData,
    integrated_adata: ad.AnnData,
    rng: np.random.Generator,
) -> ad.AnnData:
    """Create a shuffled union with aligned, intentionally missing obs metadata."""
    combined = ad.concat([base_adata, integrated_adata], join="outer")
    combined.obsm.clear()
    if "label:scANVI" in _obs(combined):
        combined.obs = _obs(combined).drop(columns="label:scANVI")
    # The modules must align by observation index, not assume source order.
    combined = combined[rng.permutation(combined.n_obs)].copy()
    return combined


def build_fixtures(
    output_dir: str | Path,
    seed: int = 42,
    *,
    writer: Callable[[Any, Path], None] | None = None,
) -> list[str]:
    """Build the six consumed H5ADs, optionally through the shared writer."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    produced: list[str] = []
    for algorithm in ("scvi", "scanvi"):
        labels = algorithm == "scanvi"
        base = create_base_adata(rng, algorithm, labels)
        integrated = create_integrated_adata(rng, labels)
        combined = create_combined_adata(base, integrated, rng)
        for kind, adata in (
            ("base", base),
            ("integrated", integrated),
            ("combined", combined),
        ):
            filename = f"{kind}_{algorithm}.h5ad"
            path = output_dir / filename
            if path.exists():
                raise FileExistsError(f"Refusing to overwrite existing fixture: {path}")
            if writer is None:
                adata.write_h5ad(path, compression="gzip", compression_opts=4)
            else:
                writer(adata, path)
            produced.append(filename)
    return produced


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    names = build_fixtures(args.output_dir, args.seed)
    print(f"Built {len(names)} fixtures in {args.output_dir}")
    for name in names:
        print(name)


if __name__ == "__main__":
    main()
