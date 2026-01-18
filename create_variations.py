#!/usr/bin/env python3
"""
Create variations of anndata files for testing the ADATA_UNIFY module.
Based on TEST_STRUCTURES.md requirements.
"""

import anndata as ad
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.sparse import csr_matrix, csc_matrix
import sys

OUTDIR = Path('anndata-variations')

OUTDIR.mkdir(parents=True, exist_ok=True)

# Load base anndata file
input_file = Path("pbmc/SRR28679756_filtered_matrix.h5ad")
if not input_file.exists():
    print(f"Error: Input file {input_file} not found!")
    sys.exit(1)

print("Loading base anndata file...")
adata_base = ad.read_h5ad(input_file)
adata_base.var["symbols"] = adata_base.var.index
adata_base.layers["counts"] = adata_base.X.copy()
print(adata_base)


# Helper function to save variation
def save_variation(adata, filename, description=""):
    """Save anndata variation to file."""
    filepath = OUTDIR / filename
    adata.write_h5ad(filepath, compression="gzip")
    if description:
        print(f"  ✓ {filename}: {description}")

# Get random cell type labels for testing
np.random.seed(42)
cell_types = ["T-Cell", "B cell", "CD4+ T", "CD8+ T", "NK cell", "Monocyte", "Dendritic", "Unassigned"]
n_cells = adata_base.n_obs

print("\n=== Creating variations ===\n")

# ============================================================================
# 1. COUNTS LAYER VARIATIONS
# ============================================================================
print("1. Counts Layer Variations")

# 1.1 Counts in X (base case - already has counts in layer, move to X)
adata = adata_base.copy()
adata.X = adata.layers["counts"].copy()
adata.layers = {}
save_variation(adata, "counts_in_X.h5ad", "Counts in X matrix")

# 1.2 Counts in layer (base case - already has this)
adata = adata_base.copy()
save_variation(adata, "counts_in_layer.h5ad", "Counts in 'counts' layer")

# 1.3 Counts in different layer name
adata = adata_base.copy()
adata.layers["raw"] = adata.layers["counts"].copy()
del adata.layers["counts"]
save_variation(adata, "counts_in_raw_layer.h5ad", "Counts in 'raw' layer")

# ============================================================================
# 2. GENE SYMBOL COLUMN VARIATIONS
# ============================================================================
print("\n2. Gene Symbol Column Variations")

# 2.1 Symbols as index (move symbols to index)
adata = adata_base.copy()
adata.var.index = adata.var["symbols"]
adata.var = adata.var.drop(columns=["symbols"])
save_variation(adata, "symbols_as_index.h5ad", "Gene symbols as var.index")

# 2.2 Symbols in var column (base case - already has this)
adata = adata_base.copy()
save_variation(adata, "symbols_in_column.h5ad", "Gene symbols in var['symbols']")

# 2.3 Symbols in different column name
adata = adata_base.copy()
adata.var["gene_symbol"] = adata.var["symbols"]
del adata.var["symbols"]
save_variation(adata, "symbols_in_gene_symbol.h5ad", "Gene symbols in var['gene_symbol']")

# ============================================================================
# 3. DUPLICATE GENE HANDLING
# ============================================================================
print("\n3. Duplicate Gene Handling")

# 3.1 No duplicates (base case - should already be unique)
adata = adata_base.copy()
save_variation(adata, "no_duplicates.h5ad", "No duplicate gene names")

# 3.2 Create duplicates for testing aggregation methods
# Duplicate some genes by copying them with the same symbol
adata = adata_base.copy()
n_dup = 50  # Number of genes to duplicate
dup_indices = np.random.choice(adata.n_vars, n_dup, replace=False)

# Create new var rows with same symbols (true duplicates)
var_dup = adata.var.iloc[dup_indices].copy()
var_dup.index = [f"{idx}_dup" for idx in var_dup.index]
new_var = pd.concat([adata.var, var_dup])

# Create new X with duplicated columns
X_new = np.hstack([adata.X.toarray(), adata.X[:, dup_indices].toarray()])
adata_dup = ad.AnnData(X=csr_matrix(X_new), obs=adata.obs.copy(), var=new_var)

# Recreate layers with correct shape (extend to match new var dimension)
adata_dup.layers = {}
for layer_name, layer_data in adata.layers.items():
    layer_array = layer_data.toarray()
    layer_new = np.hstack([layer_array, layer_array[:, dup_indices]])
    adata_dup.layers[layer_name] = csr_matrix(layer_new)

# Ensure duplicates have same symbol values for proper testing
for orig_idx in dup_indices:
    orig_gene_idx = adata.var.index[orig_idx]
    dup_idx = f"{orig_gene_idx}_dup"
    adata_dup.var.loc[dup_idx, "symbols"] = adata.var.loc[orig_gene_idx, "symbols"]

save_variation(adata_dup, "with_duplicates.h5ad", f"Contains {n_dup} duplicate genes for aggregation testing")

# ============================================================================
# 4. ISOFORM AGGREGATION
# ============================================================================
print("\n4. Isoform Aggregation")

# 4.1 With isoform suffixes (add numeric suffixes to some genes)
adata = adata_base.copy()
n_isoforms = 100
isoform_indices = np.random.choice(adata.n_vars, n_isoforms, replace=False)
for idx in isoform_indices:
    gene_name = adata.var.index[idx]
    adata.var.index.values[idx] = f"{gene_name}.1"
save_variation(adata, "with_isoforms.h5ad", f"Contains {n_isoforms} genes with numeric suffixes (e.g., GENE.1)")

# ============================================================================
# 5. BATCH COLUMN VARIATIONS
# ============================================================================
print("\n5. Batch Column Variations")

# 5.1 Batch column exists with different name
adata = adata_base.copy()
adata.obs["batch_id"] = adata.obs["sample"]
del adata.obs["sample"]
save_variation(adata, "batch_different_name.h5ad", "Batch column named 'batch_id'")

# 5.2 Batch column already named "batch"
adata = adata_base.copy()
adata.obs["batch"] = adata.obs["sample"]
del adata.obs["sample"]
save_variation(adata, "batch_correct_name.h5ad", "Batch column already named 'batch'")

# 5.3 Batch column missing
adata = adata_base.copy()
del adata.obs["sample"]
save_variation(adata, "batch_missing.h5ad", "No batch column in obs")

# 5.4 Batch column conflict case (batch exists but batch_col != "batch")
adata = adata_base.copy()
adata.obs["batch"] = "conflict_batch"
save_variation(adata, "batch_conflict.h5ad", "Both 'batch' and 'sample' exist (error case)")

# ============================================================================
# 6. LABEL COLUMN VARIATIONS
# ============================================================================
print("\n6. Label Column Variations")

# 6.1 Label column exists with different name
adata = adata_base.copy()
adata.obs["cell_type"] = np.random.choice(cell_types, n_cells)
save_variation(adata, "label_different_name.h5ad", "Label column named 'cell_type'")

# 6.2 Label column already named "label"
adata = adata_base.copy()
adata.obs["label"] = np.random.choice(cell_types, n_cells)
save_variation(adata, "label_correct_name.h5ad", "Label column already named 'label'")

# 6.3 Label column missing
adata = adata_base.copy()
save_variation(adata, "label_missing.h5ad", "No label column in obs")

# 6.4 Label column with unknown values
adata = adata_base.copy()
labels = np.random.choice(cell_types, n_cells)
labels[np.random.choice(n_cells, 100, replace=False)] = "Unassigned"
adata.obs["cell_type"] = labels
save_variation(adata, "label_with_unknown.h5ad", "Label column with 'Unassigned' values (unknown_label)")

# 6.5 Label column with NaN values
adata = adata_base.copy()
labels = np.random.choice(cell_types, n_cells)
nan_indices = np.random.choice(n_cells, 50, replace=False)
labels = labels.astype(object)
labels[nan_indices] = np.nan
adata.obs["cell_type"] = labels
save_variation(adata, "label_with_nan.h5ad", "Label column with NaN values")

# 6.6 Label column with special characters
adata = adata_base.copy()
special_labels = ["T-Cell", "B cell", "CD4+ T", "CD8+ T", "NK cell", "Monocyte"]
adata.obs["cell_type"] = np.random.choice(special_labels, n_cells)
save_variation(adata, "label_special_chars.h5ad", "Labels with spaces, dashes, special chars (T-Cell, B cell, CD4+)")

# 6.7 Label column conflict case
adata = adata_base.copy()
adata.obs["label"] = "conflict_label"
adata.obs["cell_type"] = np.random.choice(cell_types, n_cells)
save_variation(adata, "label_conflict.h5ad", "Both 'label' and 'cell_type' exist (error case)")

# ============================================================================
# 6b. CONDITION COLUMN VARIATIONS
# ============================================================================
print("\n6b. Condition Column Variations")

conditions = ["healthy", "disease", "treated"]

# 6b.1 Condition column exists with different name
adata = adata_base.copy()
adata.obs["disease_state"] = np.random.choice(conditions, n_cells)
save_variation(adata, "condition_different_name.h5ad", "Condition column named 'disease_state'")

# 6b.2 Condition column already named "condition"
adata = adata_base.copy()
adata.obs["condition"] = np.random.choice(conditions, n_cells)
save_variation(adata, "condition_correct_name.h5ad", "Condition column already named 'condition'")

# 6b.3 Condition column missing
adata = adata_base.copy()
save_variation(adata, "condition_missing.h5ad", "No condition column in obs")

# 6b.4 Condition column conflict case
adata = adata_base.copy()
adata.obs["condition"] = "conflict_condition"
adata.obs["disease_state"] = np.random.choice(conditions, n_cells)
save_variation(adata, "condition_conflict.h5ad", "Both 'condition' and 'disease_state' exist (error case)")

# 6b.5 Complete test data with label and condition (for differential expression testing)
adata = adata_base.copy()
adata.obs["label"] = np.random.choice(cell_types[:6], n_cells)  # Use cell types without "Unassigned"
adata.obs["condition"] = np.random.choice(conditions, n_cells)
save_variation(adata, "with_label_and_condition.h5ad", "Contains both label and condition columns for DE testing")

# ============================================================================
# 7. MATRIX FORMAT VARIATIONS
# ============================================================================
print("\n7. Matrix Format Variations")

# 7.1 Sparse CSR (base case - likely already CSR)
adata = adata_base.copy()
adata.X = csr_matrix(adata.X)
save_variation(adata, "matrix_csr.h5ad", "X is CSR sparse matrix")

# 7.2 Sparse CSC
adata = adata_base.copy()
adata.X = csc_matrix(adata.X)
save_variation(adata, "matrix_csc.h5ad", "X is CSC sparse matrix")

# 7.3 Different dtypes
for dtype_name, dtype in [("int32", np.int32), ("int64", np.int64), ("float64", np.float64)]:
    adata = adata_base.copy()
    adata.X = adata.X.astype(dtype)
    save_variation(adata, f"matrix_{dtype_name}.h5ad", f"X is {dtype_name}")

# ============================================================================
# 8. OBSERVATION METADATA VARIATIONS
# ============================================================================
print("\n8. Observation Metadata Variations")

# 8.1 Duplicate cell names
adata = adata_base.copy()
# Create some duplicates
dup_indices = np.random.choice(n_cells, 100, replace=False)
new_names = list(adata.obs_names)
for idx in dup_indices[:50]:
    new_names[idx] = new_names[dup_indices[50]]  # Make duplicates
adata.obs_names = new_names
save_variation(adata, "obs_duplicate_names.h5ad", "Contains duplicate cell names")

# 8.2 Sample column exists with same value (base case)
adata = adata_base.copy()
save_variation(adata, "sample_exists_same.h5ad", "Sample column exists")

# 8.3 Sample column exists with different value
adata = adata_base.copy()
adata.obs["sample"] = "different_sample_id"
save_variation(adata, "sample_exists_different.h5ad", "Sample column with different value")

# 8.4 Sample column missing
adata = adata_base.copy()
del adata.obs["sample"]
save_variation(adata, "sample_missing.h5ad", "No sample column")

# ============================================================================
# 9. AnnData STRUCTURE VARIATIONS
# ============================================================================
print("\n9. AnnData Structure Variations")

# 9.1 With obsm (embeddings)
adata = adata_base.copy()
# Add some mock embeddings
n_obs = adata.n_obs
adata.obsm["X_pca"] = np.random.randn(n_obs, 50).astype(np.float32)
adata.obsm["X_umap"] = np.random.randn(n_obs, 2).astype(np.float32)
save_variation(adata, "with_obsm.h5ad", "Contains obsm embeddings (X_pca, X_umap)")

# 9.2 With varm
adata = adata_base.copy()
n_vars = adata.n_vars
adata.varm["PCs"] = np.random.randn(n_vars, 10).astype(np.float32)
save_variation(adata, "with_varm.h5ad", "Contains varm (gene-level metadata)")

# 9.3 With uns
adata = adata_base.copy()
adata.uns["neighbors"] = {"connectivities_key": "connectivities"}
adata.uns["pca"] = {"variance_ratio": np.random.rand(50)}
save_variation(adata, "with_uns.h5ad", "Contains uns metadata (neighbors, pca)")

# 9.4 With layers (base case - already has counts)
adata = adata_base.copy()
adata.layers["spliced"] = adata.X.copy()
adata.layers["unspliced"] = adata.X.copy() * 0.5
save_variation(adata, "with_layers.h5ad", "Contains multiple layers (counts, spliced, unspliced)")

# ============================================================================
# 10. COMBINED SCENARIOS
# ============================================================================
print("\n10. Combined Scenarios")

# 10.1 Minimal structure (only X, obs, var)
adata = adata_base.copy()
adata.obsm = {}
adata.varm = {}
adata.uns = {}
adata.layers = {}
del adata.obs["sample"]
save_variation(adata, "minimal_structure.h5ad", "Minimal structure: only X, obs, var")

# 10.2 Complex structure (all components)
adata = adata_base.copy()
# Add all optional components
n_obs, n_vars = adata.shape
adata.obsm["X_pca"] = np.random.randn(n_obs, 50).astype(np.float32)
adata.obsm["X_umap"] = np.random.randn(n_obs, 2).astype(np.float32)
adata.varm["PCs"] = np.random.randn(n_vars, 10).astype(np.float32)
adata.uns["neighbors"] = {"connectivities_key": "connectivities"}
adata.uns["pca"] = {"variance_ratio": np.random.rand(50)}
adata.layers["spliced"] = adata.X.copy()
adata.layers["unspliced"] = adata.X.copy() * 0.5
# Add various obs columns
adata.obs["batch_id"] = np.random.choice(["batch1", "batch2", "batch3"], n_obs)
adata.obs["cell_type"] = np.random.choice(cell_types, n_obs)
# Add var columns
adata.var["gene_id"] = [f"ENSG{i:010d}" for i in range(n_vars)]
save_variation(adata, "complex_structure.h5ad", "Complex structure: all optional components present")

print(f"\n=== Done! Created {len(list(OUTDIR.glob('*.h5ad')))} variations in {OUTDIR}/ ===")

