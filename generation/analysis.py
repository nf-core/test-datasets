"""Produce compact, production-compatible analysis fixtures from source counts.

The source stage owns count matrices.  This stage deliberately reads only those
canonical outputs and derives the small PCA/graph inputs used by analysis tests.
It never copies historical PCA files or invents count values.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

from .core import BuildContext, ContractError, StageResult, h5ad_compression_violations, write_h5ad

SCANPY_CONTAINER = (
    "community.wave.seqera.io/library/"
    "leidenalg_python-igraph_python_pyyaml_scanpy:ce926a86f015211a"
)
SCANPY_IMAGE_ID = "sha256:2f51ddffd1da33cabd9f2b3dd8d8f0cf24e39fd70706b6926335588fa282869"
APPTAINER_PYTHON_PATH = (
    "PATH=/opt/wave/.pixi/envs/default/bin:/usr/local/sbin:/usr/local/bin:"
    "/usr/sbin:/usr/bin:/sbin:/bin"
)
APPTAINER_WORK_ENV = "MPLCONFIGDIR=/work/.mplconfig,NUMBA_CACHE_DIR=/work/.numba"
PCA_INPUTS = {
    "combined_filtered_pca.h5ad": "combined_filtered_matrix.h5ad",
    "SRR28679759_filtered_pca.h5ad": "SRR28679759_filtered_matrix.h5ad",
}


def build(context: BuildContext) -> StageResult:
    """Derive PCA, neighbours and Leiden labels from canonical source counts."""
    source_dir = context.out_dir / "pbmc"
    inputs = {output: source_dir / source for output, source in PCA_INPUTS.items()}
    missing = [str(path) for path in inputs.values() if not path.is_file()]
    if missing:
        raise ContractError("Analysis requires source-stage count outputs: " + ", ".join(missing))
    outputs = {output: source_dir / output for output in PCA_INPUTS}
    existing = [str(path) for path in outputs.values() if path.exists()]
    if existing:
        raise ContractError("Analysis refuses to overwrite PCA outputs: " + ", ".join(existing))

    observations = derive_pca_files(inputs, outputs, context.cache_dir / "analysis-work" / "pca")
    fragment_observations = _write_fragments(
        context,
        inputs["combined_filtered_pca.h5ad"],
        inputs["SRR28679759_filtered_pca.h5ad"],
    )
    produced = [path.relative_to(context.out_dir).as_posix() for path in outputs.values()]
    produced.extend(fragment_observations.pop("produced"))
    _assert_hard_limit(context, produced)
    return StageResult(
        stage="analysis",
        produced=sorted(produced),
        provenance={
            "counts_source": "canonical source-stage real count H5ADs",
            "analysis_recipe": (
                "pinned SCANPY_PCA/SCANPY_NEIGHBORS/SCANPY_LEIDEN-compatible Scanpy recipe"
            ),
            "container": SCANPY_CONTAINER,
            "container_image_id": SCANPY_IMAGE_ID,
            "pca_components": 50,
            "random_seed": 0,
            "synthetic_fragments": (
                "alignment vectors and doublet annotations are explicit structural fixtures; "
                "no scVI model or doublet detector was fabricated"
            ),
        },
        dependency_fingerprint=context.dependency_fingerprint,
        observations={"pca": observations, "fragments": fragment_observations},
    )


def derive_pca_files(
    inputs: dict[str, Path], outputs: dict[str, Path], work_root: Path | None = None
) -> dict[str, Any]:
    """Run the pinned production Scanpy image for caller-supplied count inputs.

    ``inputs`` and ``outputs`` share their output filenames.  The pinned image
    runs through Apptainer against a shared NFS staging directory.
    """
    if set(inputs) != set(outputs):
        raise ContractError("PCA input/output names must match")
    if any(path.exists() for path in outputs.values()):
        raise ContractError("Refusing to overwrite a PCA output")
    if work_root is None:
        configured = os.environ.get("SCDOWNSTREAM_CONTAINER_WORK_ROOT")
        if not configured:
            raise ContractError(
                "Set SCDOWNSTREAM_CONTAINER_WORK_ROOT for standalone PCA derivation"
            )
        work_root = Path(configured) / "analysis-pca"
    work_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="derive-", dir=work_root) as temporary:
        work = Path(temporary)
        for name, source in inputs.items():
            shutil.copy2(source, work / f"input-{name}")
        script = work / "derive.py"
        script.write_text(_DERIVATION_SCRIPT, encoding="utf-8")
        command = [
            "apptainer", "exec", "--containall", "--home", f"{work}:/work",
            "--bind", f"{work}:/work",
            "--env", APPTAINER_PYTHON_PATH,
            "--env", APPTAINER_WORK_ENV,
            "--pwd", "/work", _apptainer_image(SCANPY_CONTAINER),
            "python", "/work/derive.py", *inputs.keys(),
        ]
        try:
            subprocess.run(command, check=True, text=True, capture_output=True)
        except FileNotFoundError as error:
            raise ContractError("Apptainer is required for pinned analysis derivation") from error
        except subprocess.CalledProcessError as error:
            message = f"Pinned Scanpy derivation failed: {error.stderr.strip()}"
            raise ContractError(message) from error
        for name, destination in outputs.items():
            import anndata as ad

            # The pinned recipe writes the logical analysis object.  Publication
            # always goes through the shared writer so old readers get legacy
            # strings and identical X/raw_counts datasets are stored once.
            write_h5ad(ad.read_h5ad(work / name), destination)
    return {name: _h5ad_observation(path) for name, path in outputs.items()}


def _write_fragments(context: BuildContext, combined: Path, source59: Path) -> dict[str, Any]:
    """Create only consumer-required indexed fragment fixtures.

    These are intentionally small alignment contracts.  The embedding vectors
    are deterministic test vectors derived from a real-count subset, not a
    claimed biological scVI latent space.
    """
    output = context.out_dir / "intermediates"
    fragments = output / "fragments"
    doublets = output / "doublet_detection"
    if fragments.exists() or doublets.exists():
        raise ContractError("Analysis refuses existing fragment or doublet fixture directories")
    work_root = context.cache_dir / "analysis-work" / "fragments"
    work_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="derive-", dir=work_root) as temporary:
        work = Path(temporary)
        shutil.copy2(combined, work / "combined.h5ad")
        shutil.copy2(source59, work / "source59.h5ad")
        script = work / "fragments.py"
        script.write_text(_FRAGMENT_SCRIPT, encoding="utf-8")
        command = [
            "apptainer", "exec", "--containall", "--home", f"{work}:/work",
            "--bind", f"{work}:/work",
            "--env", APPTAINER_PYTHON_PATH,
            "--env", APPTAINER_WORK_ENV,
            "--pwd", "/work", _apptainer_image(SCANPY_CONTAINER),
            "python", "/work/fragments.py",
        ]
        try:
            subprocess.run(command, check=True, text=True, capture_output=True)
        except subprocess.CalledProcessError as error:
            message = f"Analysis fragment derivation failed: {error.stderr.strip()}"
            raise ContractError(message) from error
        fragments.mkdir(parents=True)
        doublets.mkdir(parents=True)
        for path in (work / "fragments").iterdir():
            destination = fragments / path.name
            if path.suffix == ".h5ad":
                import anndata as ad

                write_h5ad(ad.read_h5ad(path), destination)
            else:
                shutil.copy2(path, destination)
        for path in (work / "doublets").iterdir():
            shutil.copy2(path, doublets / path.name)
    produced = [
        path.relative_to(context.out_dir).as_posix()
        for directory in (fragments, doublets)
        for path in directory.rglob("*")
        if path.is_file()
    ]
    return {
        "produced": sorted(produced),
        "base_shape": [16, 32],
        "partial_rows": 12,
        "doublet_rows": _obs_count(source59),
    }


def _h5ad_observation(path: Path) -> dict[str, Any]:
    import anndata as ad

    adata = ad.read_h5ad(path, backed="r")
    try:
        required_layers = {"raw_counts"}
        missing = required_layers - set(adata.layers.keys())
        if missing or "X_pca" not in adata.obsm or "leiden" not in adata.obs:
            raise ContractError(f"Derived analysis H5AD is incomplete: {path}")
        if adata.obsm["X_pca"].shape[1] != 50:
            raise ContractError(f"Derived PCA does not have 50 components: {path}")
        violations = h5ad_compression_violations(path)
        if violations:
            message = f"Derived H5AD lacks gzip compression: {path}: {', '.join(violations)}"
            raise ContractError(message)
        return {
            "shape": list(adata.shape),
            "bytes": path.stat().st_size,
            "gzip": True,
            "layers": sorted(adata.layers.keys()),
        }
    finally:
        adata.file.close()


def _obs_count(path: Path) -> int:
    import anndata as ad

    adata = ad.read_h5ad(path, backed="r")
    try:
        return int(adata.n_obs)
    finally:
        adata.file.close()


def _assert_hard_limit(context: BuildContext, produced: list[str]) -> None:
    hard_limit = context.contracts["size_policy"]["hard_limit_bytes"]
    too_large = [
        relative
        for relative in produced
        if (context.out_dir / relative).stat().st_size >= hard_limit
    ]
    if too_large:
        raise ContractError("Analysis stage exceeded hard size limit: " + ", ".join(too_large))


def _apptainer_image(image: str) -> str:
    """Use an explicit registry transport instead of Docker runtime execution."""
    return image if "://" in image else f"docker://{image}"


_DERIVATION_SCRIPT = r'''import sys
import scanpy as sc
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

sc.settings.n_jobs = 1
def legacy_strings(adata):
    for frame in (adata.obs, adata.var):
        frame.index = pd.Index(frame.index.to_numpy(dtype=object), dtype=object)
        for column in frame.columns:
            if isinstance(frame[column].dtype, pd.StringDtype):
                frame[column] = frame[column].astype(object).where(frame[column].notna(), np.nan)

with threadpool_limits(limits=1):
    for name in sys.argv[1:]:
        adata = sc.read_h5ad("input-" + name)
        # X remains source counts. ``raw_counts`` is the sole retained alias:
        # it is the canonical count-provenance contract. ``counts`` and the
        # temporary PCA normalisation layer are not fixture-consumer inputs.
        source_counts = adata.X.copy()
        for layer in list(adata.layers):
            del adata.layers[layer]
        adata.X = source_counts
        adata.layers["raw_counts"] = source_counts.copy()
        adata.layers["pca_log1p"] = source_counts.copy()
        sc.pp.normalize_total(adata, layer="pca_log1p")
        sc.pp.log1p(adata, layer="pca_log1p")
        # Keep the full >2,000 real-gene panel available to feature-selection
        # branches; PCA itself follows the production SCANPY_PCA recipe.
        adata.var["highly_variable"] = True
        sc.pp.pca(adata, n_comps=50, layer="pca_log1p", random_state=0, key_added="X_pca")
        # PCA loadings are ~600 KB dense float64 matrices.  The canonical PCA
        # consumers use scores, not loadings; PCA/module tests create and test
        # their own loadings from count-only source fixtures.
        adata.varm.clear()
        sc.pp.neighbors(
            adata,
            n_neighbors=15,
            use_rep="X_pca",
            method="umap",
            metric="euclidean",
            random_state=0,
        )
        adata.obsp["distances"] = adata.obsp["distances"].astype("float32")
        sc.tl.leiden(adata, resolution=1.0, random_state=0, flavor="igraph", key_added="leiden")
        # Scanpy records ``base=None`` after log1p, which AnnData serialises as
        # an HDF5 null dataset.  The shared gzip contract intentionally has no
        # null-dataset exemption; the layer itself remains self-describing.
        adata.uns.pop("log1p", None)
        del adata.layers["pca_log1p"]
        # Individual source matrices intentionally avoid canonical design
        # columns before ADATA_UNIFY.  PCA-only consumers do not pass through
        # that stage, so expose read-only aliases only in this derivative.
        for canonical, synthetic in (
            ("donor", "synthetic_donor"),
            ("condition", "synthetic_condition"),
            ("label", "synthetic_label"),
        ):
            if canonical not in adata.obs and synthetic in adata.obs:
                adata.obs[canonical] = adata.obs[synthetic].copy()
        legacy_strings(adata)
        adata.write_h5ad(name, compression="gzip", compression_opts=4)
'''


_FRAGMENT_SCRIPT = r'''import pickle
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd

# The pipeline still has pandas readers predating StringDtype pickle support.
# Keep these interchange tables on the stable object-string representation.
pd.options.future.infer_string = False

fragments = Path("fragments")
doublets = Path("doublets")
fragments.mkdir()
doublets.mkdir()
base = ad.read_h5ad("combined.h5ad")[:16, :32].copy()
for frame in (base.obs, base.var):
    frame.index = pd.Index(frame.index.to_numpy(dtype=object), dtype=object)
    for column in frame.columns:
        if isinstance(frame[column].dtype, pd.StringDtype):
            frame[column] = frame[column].astype(object).where(frame[column].notna(), np.nan)
base.write_h5ad(fragments / "merged_outer.h5ad", compression="gzip", compression_opts=4)
# Pickles cross the pinned process environments.  Avoid pandas' versioned
# StringDtype serialisation by using ordinary object barcode indexes.
index = pd.Index(base.obs_names.astype(str).to_numpy(dtype=object), dtype=object)
# Explicit deterministic structural vectors, calculated from the retained real
# count subset only to make alignment easy to inspect.  They are not scVI output.
counts = base.X.toarray() if hasattr(base.X, "toarray") else np.asarray(base.X)
weights = np.arange(1, 9, dtype=np.float32)[None, :]
vectors = (counts[:, :8].astype(np.float32) + weights) / (
    1.0 + counts[:, :8].sum(axis=1, keepdims=True)
)
pd.DataFrame(vectors, index=index).to_pickle(fragments / "X_scvi.pkl")
pd.DataFrame(vectors[:12, :2], index=index[:12]).to_pickle(fragments / "X_scvi-global_umap.pkl")
pd.DataFrame(
    {"scvi-global-1.0_leiden": pd.Categorical(["0", "1"] * 6)}, index=index[:12]
).to_pickle(fragments / "scvi-global-1.0_leiden.pkl")
pd.DataFrame(
    {"scvi-global-1.0:entropy": np.linspace(0.0, 1.0, 12, dtype=np.float32)}, index=index[:12]
).to_pickle(fragments / "scvi-global-1.0_entropy.pkl")
# ADATA_EXTEND consumes this pickle as an uns value; it does not require a
# graph algorithm result for this structural merge test.
with open(fragments / "scvi-global-1.0_paga.pkl", "wb") as handle:
    pickle.dump(
        {"groups": "scvi-global-1.0_leiden", "params": {"groups": "scvi-global-1.0_leiden"}},
        handle,
    )

source = ad.read_h5ad("source59.h5ad")
source_index = pd.Index(source.obs_names.astype(str).to_numpy(dtype=object), dtype=object)
for name, modulus, offset in (
    ("test_doubletdetection", 7, 0), ("test_scrublet", 9, 0), ("test_solo", 11, 0)
):
    values = np.fromiter(
        ((position + offset) % modulus == 0 for position in range(source.n_obs)), dtype=bool
    )
    pd.DataFrame({name: values}, index=source_index).to_pickle(doublets / f"{name}.pkl")
'''
