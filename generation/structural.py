"""Structural, embedding-alignment and plotting fixture producer."""

from __future__ import annotations

import importlib.util
import time
from pathlib import Path
from typing import Any

from .core import BuildContext, StageResult, write_h5ad


def _load_recipe(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load fixture recipe: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build(context: BuildContext) -> StageResult:
    """Build only the structural stage's declared files using shared H5AD writing."""
    started = time.perf_counter()
    repository = context.contracts_path.parent
    variations = _load_recipe(
        "scdownstream_create_variations", repository / "create_variations.py"
    )
    embeddings = _load_recipe(
        "scdownstream_merge_embedding_fixtures",
        repository / "intermediates/merge_embeddings/generate.py",
    )
    volcano = _load_recipe(
        "scdownstream_volcano_fixtures", repository / "volcanoplot/make_fixtures.py"
    )

    produced: list[str] = []
    observations: dict[str, Any] = {"fixture_details": {}}
    structural_dir = context.out_dir / "anndata-variations"
    structural_manifest = variations.build_variations(
        structural_dir,
        seed=context.seed,
        writer=write_h5ad,
        write_manifest=False,
    )
    for filename, details in structural_manifest["files"].items():
        relative = f"anndata-variations/{filename}"
        produced.append(relative)
        observations["fixture_details"][relative] = details

    embedding_dir = context.out_dir / "intermediates/merge_embeddings"
    for filename in embeddings.build_fixtures(
        embedding_dir, seed=context.seed, writer=write_h5ad
    ):
        relative = f"intermediates/merge_embeddings/{filename}"
        produced.append(relative)
    observations["embedding_alignment"] = {
        "files": embeddings.FILENAMES,
        "combined_obs": 80,
        "combined_order": "seeded permutation of the 50 base and 30 integrated observations",
        "metadata_alignment": (
            "base_only and integrated_only are source-specific and NaN outside their source"
        ),
        "algorithms": ["scvi", "scanvi"],
    }

    volcano_dir = context.out_dir / "volcanoplot"
    for filename in volcano.build_fixtures(volcano_dir, seed=context.seed):
        relative = f"volcanoplot/{filename}"
        produced.append(relative)
    observations["volcanoplot_files"] = volcano.FILENAMES
    observations["synthetic_counts"] = True
    observations["raw_counts_layer_written"] = False
    observations["seed"] = context.seed

    return StageResult(
        stage="structural",
        produced=produced,
        provenance={
            "recipes": [
                "create_variations.py",
                "intermediates/merge_embeddings/generate.py",
                "volcanoplot/make_fixtures.py",
            ],
            "data_source": "synthetic structural/test-edge-case data; no biological source files",
        },
        dependency_fingerprint=context.dependency_fingerprint,
        observations=observations,
        elapsed_seconds=time.perf_counter() - started,
    )
