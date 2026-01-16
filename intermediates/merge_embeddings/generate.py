#!/usr/bin/env python3
"""
Generate test data for adata/mergeembeddings module.

This script creates three h5ad files that simulate the merge embeddings scenario:
- base.h5ad: Pre-existing data with X_{integration} embeddings
- integrated.h5ad: Newly integrated data with X_emb embeddings
- combined.h5ad: Union of cells from both (no embeddings yet)

For scanvi tests, additional files include label:scANVI columns.
"""

import numpy as np
import anndata as ad
import pandas as pd
from scipy.sparse import csr_matrix
from pathlib import Path

np.random.seed(42)

OUTPUT_DIR = Path(__file__).parent
N_BASE_CELLS = 50
N_INTEGRATED_CELLS = 30
N_GENES = 100
EMBEDDING_DIM = 10


def create_base_adata(integration_name: str = "scvi", include_scanvi_labels: bool = False):
    """Create base anndata with X_{integration} embedding."""
    cell_names = [f"base_cell_{i}" for i in range(N_BASE_CELLS)]
    gene_names = [f"gene_{i}" for i in range(N_GENES)]

    # Random count matrix
    X = csr_matrix(np.random.poisson(5, size=(N_BASE_CELLS, N_GENES)).astype(np.float32))

    # Create embedding
    embedding = np.random.randn(N_BASE_CELLS, EMBEDDING_DIM).astype(np.float32)

    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=cell_names),
        var=pd.DataFrame(index=gene_names),
    )
    adata.obsm[f"X_{integration_name}"] = embedding

    if include_scanvi_labels:
        labels = np.random.choice(["TypeA", "TypeB", "TypeC"], size=N_BASE_CELLS)
        adata.obs["label:scANVI"] = pd.Categorical(labels)

    return adata


def create_integrated_adata(include_scanvi_labels: bool = False):
    """Create integrated anndata with X_emb embedding."""
    cell_names = [f"integrated_cell_{i}" for i in range(N_INTEGRATED_CELLS)]
    gene_names = [f"gene_{i}" for i in range(N_GENES)]

    # Random count matrix
    X = csr_matrix(np.random.poisson(5, size=(N_INTEGRATED_CELLS, N_GENES)).astype(np.float32))

    # Create embedding (named X_emb as expected by the module)
    embedding = np.random.randn(N_INTEGRATED_CELLS, EMBEDDING_DIM).astype(np.float32)

    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=cell_names),
        var=pd.DataFrame(index=gene_names),
    )
    adata.obsm["X_emb"] = embedding

    if include_scanvi_labels:
        labels = np.random.choice(["TypeA", "TypeB", "TypeC"], size=N_INTEGRATED_CELLS)
        adata.obs["label:scANVI"] = pd.Categorical(labels)

    return adata


def create_combined_adata(base_adata: ad.AnnData, integrated_adata: ad.AnnData):
    """Create combined anndata (union of cells, no embeddings)."""
    # Concatenate the two datasets
    combined = ad.concat([base_adata, integrated_adata], join="outer")

    # Remove embeddings (the module will add them)
    combined.obsm = {}

    # Remove scanvi labels if present (the module will add them for scanvi)
    if "label:scANVI" in combined.obs.columns:
        del combined.obs["label:scANVI"]

    return combined


def main():
    print(f"Generating test data in {OUTPUT_DIR}")

    # Standard scvi test data
    print("Creating scvi test data...")
    base_scvi = create_base_adata(integration_name="scvi", include_scanvi_labels=False)
    integrated_scvi = create_integrated_adata(include_scanvi_labels=False)
    combined_scvi = create_combined_adata(base_scvi, integrated_scvi)

    base_scvi.write_h5ad(OUTPUT_DIR / "base_scvi.h5ad")
    integrated_scvi.write_h5ad(OUTPUT_DIR / "integrated_scvi.h5ad")
    combined_scvi.write_h5ad(OUTPUT_DIR / "combined_scvi.h5ad")

    # scanvi test data (includes label:scANVI)
    print("Creating scanvi test data...")
    base_scanvi = create_base_adata(integration_name="scanvi", include_scanvi_labels=True)
    integrated_scanvi = create_integrated_adata(include_scanvi_labels=True)
    combined_scanvi = create_combined_adata(base_scanvi, integrated_scanvi)

    base_scanvi.write_h5ad(OUTPUT_DIR / "base_scanvi.h5ad")
    integrated_scanvi.write_h5ad(OUTPUT_DIR / "integrated_scanvi.h5ad")
    combined_scanvi.write_h5ad(OUTPUT_DIR / "combined_scanvi.h5ad")

    print("Done! Created files:")
    for f in sorted(OUTPUT_DIR.glob("*.h5ad")):
        print(f"  - {f.name}")


if __name__ == "__main__":
    main()
