"""Shared contracts, serialisation and stage lifecycle support."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from time import perf_counter
from typing import Any, Literal

StageName = Literal["source", "structural", "analysis", "reference", "singler"]
STAGES: tuple[StageName, ...] = ("source", "structural", "analysis", "reference", "singler")


class ContractError(ValueError):
    """Raised when an authored contract cannot safely drive a build."""


class StageNotImplementedError(NotImplementedError):
    """Raised by a declared producer before it can produce valid artefacts."""


@dataclass(frozen=True)
class BuildContext:
    """Inputs shared by every producer stage."""

    stage: StageName
    contracts: dict[str, Any]
    contracts_path: Path
    out_dir: Path
    cache_dir: Path
    pipeline_dir: Path
    seed: int
    pipeline_revision: str
    dependency_fingerprint: str
    started_at: str


@dataclass(frozen=True)
class StageResult:
    """A producer's publication-ready account of its output."""

    stage: StageName
    produced: list[str]
    provenance: dict[str, Any]
    dependency_fingerprint: str
    observations: dict[str, Any] = field(default_factory=dict)
    elapsed_seconds: float = 0.0

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def load_contracts(path: Path) -> dict[str, Any]:
    """Load and validate the small authored contract surface used by all stages."""
    try:
        contracts = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as error:
        raise ContractError(f"Contract file does not exist: {path}") from error
    except json.JSONDecodeError as error:
        raise ContractError(f"Contract file is not valid JSON: {path}: {error}") from error

    required = {"schema_version", "seed", "pipeline", "size_policy", "stages", "fixture_families"}
    missing = sorted(required - contracts.keys())
    if missing:
        raise ContractError(f"Contract file is missing keys: {', '.join(missing)}")
    if contracts["schema_version"] != 1:
        raise ContractError(f"Unsupported contract schema: {contracts['schema_version']}")
    if not isinstance(contracts["seed"], int):
        raise ContractError("Contract seed must be an integer")
    if set(contracts["stages"]) != set(STAGES):
        raise ContractError(
            "Contract stages must define source, structural, analysis, reference and singler"
        )
    size_policy = contracts["size_policy"]
    if size_policy.get("hard_limit_bytes", 0) <= size_policy.get("target_bytes", 0):
        raise ContractError("The hard size limit must exceed the target size")
    return contracts


def sha256_file(path: Path) -> str:
    """Return the SHA-256 checksum of a file without loading it into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pipeline_revision(pipeline_dir: Path) -> str:
    """Resolve the checked-out pipeline revision required by the contract."""
    try:
        return subprocess.check_output(
            ["git", "-C", str(pipeline_dir), "rev-parse", "HEAD"], text=True
        ).strip()
    except subprocess.CalledProcessError as error:
        raise ContractError(f"Cannot resolve pipeline revision in {pipeline_dir}") from error


def dependency_fingerprint(contracts_path: Path, pipeline_dir: Path, revision: str) -> str:
    """Fingerprint immutable inputs that make a cached stage result reusable."""
    repository = contracts_path.parent
    recipes = [repository / "requirements.txt"]
    recipes.extend(sorted((repository / "generation").glob("*.py")))
    recipes.extend(sorted((repository / "generation").glob("*.R")))
    recipes.append(repository / "create_variations.py")
    recipes.append(repository / "pbmc/generate_pca_fixtures.py")
    recipes.append(repository / "intermediates/merge_embeddings/generate.py")
    recipes.append(repository / "volcanoplot/make_fixtures.py")
    tracked_paths = subprocess.check_output(
        ["git", "-C", str(pipeline_dir), "ls-files", "-z"], text=True
    ).split("\0")
    pipeline_files = {
        relative: sha256_file(pipeline_dir / relative)
        for relative in tracked_paths
        if relative
        and (pipeline_dir / relative).is_file()
        and not relative.endswith((".snap", ".png", ".pdf"))
    }
    payload = {
        "contracts_sha256": sha256_file(contracts_path),
        "pipeline_dir": str(pipeline_dir.resolve()),
        "pipeline_revision": revision,
        "pipeline_files": pipeline_files,
        "python": sys.version,
        "recipes": {
            path.relative_to(repository).as_posix(): sha256_file(path)
            for path in recipes
            if path.is_file()
        },
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def create_context(
    *,
    stage: StageName,
    contracts_path: Path,
    out_dir: Path,
    cache_dir: Path,
    pipeline_dir: Path,
) -> BuildContext:
    """Resolve a safe build context without creating or replacing output files."""
    contracts_path = contracts_path.resolve()
    contracts = load_contracts(contracts_path)
    pipeline_dir = pipeline_dir.resolve()
    out_dir = out_dir.resolve()
    cache_dir = cache_dir.resolve()
    repository_dir = contracts_path.parent
    if _is_relative_to(out_dir, repository_dir) or _is_relative_to(out_dir, pipeline_dir):
        raise ContractError("--out must be outside the dataset and pipeline checkouts")
    revision = pipeline_revision(pipeline_dir)
    expected_revision = contracts["pipeline"]["revision"]
    if revision != expected_revision:
        raise ContractError(
            f"Pipeline revision is {revision}, but contracts require {expected_revision}"
        )
    return BuildContext(
        stage=stage,
        contracts=contracts,
        contracts_path=contracts_path,
        out_dir=out_dir,
        cache_dir=cache_dir,
        pipeline_dir=pipeline_dir,
        seed=contracts["seed"],
        pipeline_revision=revision,
        dependency_fingerprint=dependency_fingerprint(contracts_path, pipeline_dir, revision),
        started_at=datetime.now(UTC).isoformat(),
    )


def stage_not_implemented(stage: str, context: BuildContext) -> StageNotImplementedError:
    """Describe an unavailable stage without ever representing it as successful."""
    families = ", ".join(context.contracts["stages"][stage])
    return StageNotImplementedError(
        f"Stage '{stage}' is not implemented. It owns: {families}. No output was written."
    )


def assert_stage_output_is_available(context: BuildContext) -> None:
    """Reject collisions before a producer can overwrite another stage's result."""
    stage_families = context.contracts["stages"]
    owned = [
        context.out_dir / family
        for family in stage_families[context.stage]
        if sum(family in families for families in stage_families.values()) == 1
    ]
    owned.append(context.out_dir / "reports" / f"{context.stage}.json")
    collisions = [str(path) for path in owned if path.exists()]
    if collisions:
        raise ContractError(
            f"Stage '{context.stage}' refuses existing owned output: {', '.join(collisions)}"
        )


def write_stage_report(context: BuildContext, result: StageResult, started: float) -> Path:
    """Write a producer result only after the producer has completed successfully."""
    if result.stage != context.stage:
        raise ContractError(f"Stage '{context.stage}' returned a result for '{result.stage}'")
    if result.dependency_fingerprint != context.dependency_fingerprint:
        raise ContractError("Stage result dependency fingerprint does not match its build context")
    if not result.produced or len(set(result.produced)) != len(result.produced):
        raise ContractError("A successful stage must report a non-empty, unique output list")
    hard_limit = context.contracts["size_policy"]["hard_limit_bytes"]
    previously_produced: set[str] = set()
    for prior_report in (context.out_dir / "reports").glob("*.json"):
        prior = json.loads(prior_report.read_text(encoding="utf-8"))
        previously_produced.update(prior.get("produced", []))
    for relative in result.produced:
        path = (context.out_dir / relative).resolve()
        if not _is_relative_to(path, context.out_dir):
            raise ContractError(f"Stage output escapes candidate directory: {relative}")
        if relative in previously_produced:
            raise ContractError(f"Stage output already belongs to another stage: {relative}")
        if not path.is_file():
            raise ContractError(f"Stage output is missing: {relative}")
        if path.stat().st_size >= hard_limit:
            raise ContractError(f"Stage output exceeds the hard size limit: {relative}")
        if path.suffix == ".h5ad" and h5ad_compression_violations(path):
            raise ContractError(f"Stage output has uncompressed H5AD arrays: {relative}")
    report = context.out_dir / "reports" / f"{context.stage}.json"
    if report.exists():
        raise ContractError(f"Stage report already exists: {report}")
    report.parent.mkdir(parents=True, exist_ok=True)
    payload = result.to_dict()
    payload["elapsed_seconds"] = perf_counter() - started
    payload["started_at"] = context.started_at
    payload["pipeline_revision"] = context.pipeline_revision
    payload["artefacts"] = {
        relative: {
            "sha256": sha256_file(context.out_dir / relative),
            "bytes": (context.out_dir / relative).stat().st_size,
        }
        for relative in result.produced
    }
    report.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def write_h5ad(adata: Any, path: Path) -> None:
    """Write an AnnData object with the required internal gzip representation."""
    import h5py
    import numpy as np
    import pandas as pd

    if path.exists():
        raise ContractError(f"Refusing to overwrite an existing H5AD: {path}")
    adata = adata.copy()
    # Older production readers cannot write pandas nullable-string encodings.
    # Preserve missing values; do not turn them into literal 'None' or '<NA>'.
    for frame in (adata.obs, adata.var):
        frame.index = pd.Index(frame.index.to_numpy(dtype=object), dtype=object)
        for column in frame.columns:
            if isinstance(frame[column].dtype, pd.StringDtype):
                frame[column] = frame[column].astype(object).where(frame[column].notna(), np.nan)
    path.parent.mkdir(parents=True, exist_ok=True)
    shared_counts = False
    if (
        "raw_counts" in adata.layers
        and adata.X is not None
        and type(adata.X) is type(adata.layers["raw_counts"])
        and adata.X.dtype == adata.layers["raw_counts"].dtype
    ):
        difference = adata.X != adata.layers["raw_counts"]
        unequal = difference.nnz if hasattr(difference, "nnz") else np.any(difference)
        if not unequal:
            # Omit the duplicate before writing: deleting it afterwards leaves
            # allocated free space in the HDF5 file instead of shrinking it.
            del adata.layers["raw_counts"]
            shared_counts = True
    adata.write_h5ad(path, compression="gzip", compression_opts=4)
    if shared_counts:
        with h5py.File(path, "r+") as handle:
            handle["layers"]["raw_counts"] = handle["X"]
    violations = h5ad_compression_violations(path)
    if violations:
        path.unlink(missing_ok=True)
        raise ContractError(
            f"H5AD writer did not gzip all arrays in {path}: {', '.join(violations)}"
        )


def rewrite_h5ad_with_gzip(source: Path, destination: Path) -> None:
    """Losslessly rewrite an H5AD through the common writer."""
    import anndata as ad

    write_h5ad(ad.read_h5ad(source), destination)


def h5ad_compression_violations(path: Path) -> list[str]:
    """Return non-empty, non-scalar HDF5 datasets lacking gzip compression."""
    import h5py

    violations: list[str] = []
    with h5py.File(path, "r") as handle:

        def inspect(name: str, node: Any) -> None:
            if not isinstance(node, h5py.Dataset):
                return
            if node.shape is None or node.shape == () or 0 in node.shape:
                return
            if node.compression != "gzip":
                violations.append(name)

        handle.visititems(inspect)
    return violations


def _is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True
