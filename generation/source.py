"""Build count-preserving primary scdownstream source representations."""

from __future__ import annotations

import hashlib
import shutil
import subprocess
import tempfile
import urllib.request
from pathlib import Path
from typing import Any, cast

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from .core import BuildContext, ContractError, StageResult, h5ad_compression_violations, write_h5ad

SAMPLES = ("SRR28679756", "SRR28679759")
RDS_SAMPLES = (
    ("SRR28679757", "sce"),
    ("SRR28679758", "seurat"),
)
REQUIRED_GENES = ("ANXA1", "JCHAIN", "GPNMB", "HSPA12B", "VCAN")
RAW_BACKGROUND_CELLS = 512
RAW_SHALLOW_BACKGROUND_CELLS = 256
CSV_CELLS = 64
CSV_GENES = 128


def build(context: BuildContext) -> StageResult:
    """Build source formats from checksum-verified immutable source blobs."""
    source_contract = context.contracts["sources"]["pbmc"]
    sources = _verify_sources(context.cache_dir, source_contract["artifacts"])
    selected = _build_h5ad_sources(context, sources)
    _write_csv(
        selected["SRR28679756"], context.out_dir / "pbmc" / "SRR28679756_filtered_matrix.csv"
    )
    _write_10x_h5(
        selected["SRR28679759"], context.out_dir / "misc" / "filtered_feature_bc_matrix.h5"
    )
    rds_observations = _build_rds_sources(context, sources)
    reader_observations = _verify_production_readers(context)
    produced = sorted(
        path.relative_to(context.out_dir).as_posix()
        for root in (context.out_dir / "pbmc", context.out_dir / "misc")
        for path in root.rglob("*")
        if path.is_file()
    )
    _enforce_size_limit(context, produced)
    return StageResult(
        stage="source",
        produced=produced,
        provenance={
            "source_revision": source_contract["revision"],
            "source_artifacts": {
                name: {"sha256": spec["sha256"], "url": spec["url"]}
                for name, spec in source_contract["artifacts"].items()
            },
            "ena_accessions": source_contract["ena_accessions"],
            "processing_history": source_contract["processing_history"],
            "runtime": source_contract["r_runtime"],
        },
        dependency_fingerprint=context.dependency_fingerprint,
        observations={
            "selection": source_contract["selection"],
            "h5ad": _h5ad_observations(context.out_dir, produced),
            "rds": rds_observations,
            "production_readers": reader_observations,
        },
    )


def _verify_sources(cache_dir: Path, artifacts: dict[str, dict[str, str]]) -> dict[str, Path]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    verified: dict[str, Path] = {}
    for name, spec in artifacts.items():
        path = cache_dir / name
        if path.exists() and _sha256(path) != spec["sha256"]:
            raise ContractError(f"Cached source checksum does not match contract: {path}")
        if not path.exists():
            temporary = path.with_suffix(path.suffix + ".partial")
            try:
                with (
                    urllib.request.urlopen(spec["url"]) as response,
                    temporary.open("wb") as output,
                ):
                    shutil.copyfileobj(response, output)
                if _sha256(temporary) != spec["sha256"]:
                    raise ContractError(
                        f"Downloaded source checksum does not match contract: {name}"
                    )
                temporary.replace(path)
            finally:
                temporary.unlink(missing_ok=True)
        verified[name] = path
    return verified


def _build_h5ad_sources(context: BuildContext, sources: dict[str, Path]) -> dict[str, ad.AnnData]:
    filtered = {
        sample: ad.read_h5ad(sources[f"{sample}_filtered_matrix.h5ad"]) for sample in SAMPLES
    }
    raw = {sample: ad.read_h5ad(sources[f"{sample}_raw_matrix.h5ad"]) for sample in SAMPLES}
    selected_indices = {
        sample: _deepest_indices(adata.X, adata.obs_names, 384)
        for sample, adata in filtered.items()
    }
    common_genes = _common_gene_indices(filtered, selected_indices, 3000)
    selected: dict[str, ad.AnnData] = {}
    for sample in SAMPLES:
        filtered_adata = filtered[sample]
        raw_adata = raw[sample]
        cell_indices = selected_indices[sample]
        filtered_output = filtered_adata[cell_indices, common_genes].copy()
        _annotate_filtered(filtered_output, sample)
        raw_indices, selected_mask = _matched_raw_indices(
            raw_adata, filtered_adata, cell_indices, common_genes
        )
        raw_output = raw_adata[raw_indices, common_genes].copy()
        filtered_output.var["symbols"] = filtered_output.var_names.to_numpy(dtype=object)
        raw_output.var["symbols"] = raw_output.var_names.to_numpy(dtype=object)
        _annotate_raw(raw_output, sample, filtered_output, selected_mask)
        _assert_pair_contract(filtered_output, raw_output)
        _assert_individual_source_metadata(filtered_output, sample)
        _assert_individual_source_metadata(raw_output, sample)
        write_h5ad(filtered_output, context.out_dir / "pbmc" / f"{sample}_filtered_matrix.h5ad")
        write_h5ad(raw_output, context.out_dir / "pbmc" / f"{sample}_raw_matrix.h5ad")
        selected[sample] = filtered_output

    combined = ad.concat([selected[sample] for sample in SAMPLES], join="inner", merge="same")
    _materialise_combined_design(combined)
    combined.var["highly_variable"] = True
    write_h5ad(combined, context.out_dir / "pbmc" / "combined_filtered_matrix.h5ad")
    return selected


def _deepest_indices(matrix: Any, names: pd.Index, count: int) -> np.ndarray:
    depth = _row_sums(matrix)
    if len(depth) < count:
        raise ContractError(f"Source contains {len(depth)} cells, but selection requires {count}")
    return np.lexsort((names.astype(str).to_numpy(), -depth))[:count]


def _common_gene_indices(
    filtered: dict[str, ad.AnnData], selected_indices: dict[str, np.ndarray], count: int
) -> np.ndarray:
    shared = filtered[SAMPLES[0]].var_names.astype(str)
    for sample in SAMPLES[1:]:
        if not shared.equals(filtered[sample].var_names.astype(str)):
            raise ContractError("Pinned filtered sources do not have an identical feature axis")
    selected = sparse.vstack(
        [cast(Any, filtered[sample].X)[selected_indices[sample]] for sample in SAMPLES]
    )
    scores = _variance(selected)
    names = shared.to_numpy()
    positions = {name: index for index, name in enumerate(names)}
    missing = [gene for gene in REQUIRED_GENES if gene not in positions]
    if missing:
        raise ContractError(f"Required source genes are absent: {', '.join(missing)}")
    ranked = np.lexsort((names, -scores))
    chosen: list[int] = []
    for index in [positions[gene] for gene in REQUIRED_GENES] + ranked.tolist():
        if index not in chosen:
            chosen.append(index)
        if len(chosen) == count:
            break
    return np.asarray(sorted(chosen), dtype=int)


def _matched_raw_indices(
    raw: ad.AnnData, filtered: ad.AnnData, filtered_indices: np.ndarray, gene_indices: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    raw_names = raw.obs_names.astype(str)
    selected_names = filtered.obs_names[filtered_indices].astype(str)
    selected_positions = raw_names.get_indexer(selected_names)
    if (selected_positions < 0).any():
        raise ContractError("Pinned raw source lacks selected filtered barcodes")
    # SCAR's default ambient-profile iterations need genuine near-empty
    # droplets.  Depth is evaluated on the retained feature panel, not the
    # source's discarded genes, so this remains a count-preserving subset.
    raw_depth = _row_sums(cast(Any, raw.X)[:, gene_indices])
    filtered_names = set(filtered.obs_names.astype(str))
    background = np.asarray(
        [
            index
            for index, name in enumerate(raw_names)
            if name not in filtered_names and 0 < raw_depth[index] < 100
        ],
        dtype=int,
    )
    if len(background) < RAW_BACKGROUND_CELLS:
        raise ContractError(
            "Pinned raw source lacks 512 non-zero background droplets below 100 counts"
        )
    shallow = background[raw_depth[background] == 2]
    if len(shallow) < RAW_SHALLOW_BACKGROUND_CELLS:
        raise ContractError(
            "Pinned raw source lacks 256 two-count droplets on the selected gene panel"
        )
    shallow = shallow[np.lexsort((raw_names[shallow].to_numpy(), raw_depth[shallow]))][
        :RAW_SHALLOW_BACKGROUND_CELLS
    ]
    remaining = background[~np.isin(background, shallow)]
    deep_count = RAW_BACKGROUND_CELLS - RAW_SHALLOW_BACKGROUND_CELLS
    deep = remaining[
        np.lexsort((raw_names[remaining].to_numpy(), -raw_depth[remaining]))[:deep_count]
    ]
    if len(deep) != deep_count:
        raise ContractError(
            "Pinned raw source lacks enough background droplets after shallow selection"
        )
    background = np.concatenate((shallow, deep))
    indices = np.concatenate((selected_positions, background))
    return indices, np.arange(len(indices)) < len(selected_positions)


def _annotate_filtered(adata: ad.AnnData, sample: str) -> None:
    original_names = adata.obs_names.astype(str)
    adata.obs_names = [f"{sample}::{name}" for name in original_names]
    _add_identity_and_design(adata, sample, np.ones(adata.n_obs, dtype=bool))


def _annotate_raw(
    adata: ad.AnnData, sample: str, filtered: ad.AnnData, selected_mask: np.ndarray
) -> None:
    original_names = adata.obs_names.astype(str)
    adata.obs_names = [f"{sample}::{name}" for name in original_names]
    _add_identity_and_design(adata, sample, selected_mask)
    filtered_obs = cast(pd.DataFrame, filtered.obs)
    raw_obs = cast(pd.DataFrame, adata.obs)
    selected_design = filtered_obs.set_index("source_barcode")[
        ["synthetic_donor", "synthetic_condition", "synthetic_label"]
    ]
    selected_rows = raw_obs["is_selected_filtered_cell"]
    for column in ("synthetic_donor", "synthetic_condition", "synthetic_label"):
        raw_obs.loc[selected_rows, column] = raw_obs.loc[selected_rows, "source_barcode"].map(
            selected_design[column]
        )


def _add_identity_and_design(adata: ad.AnnData, sample: str, selected_mask: np.ndarray) -> None:
    barcodes = adata.obs_names.astype(str).str.split("::", n=1).str[-1]
    obs = cast(pd.DataFrame, adata.obs)
    obs.drop(columns=["batch", "donor", "condition", "label"], errors="ignore", inplace=True)
    obs["source_sample"] = sample
    obs["source_barcode"] = barcodes.to_numpy()
    obs["sample"] = sample
    obs["is_selected_filtered_cell"] = selected_mask
    obs["source_role"] = np.where(selected_mask, "filtered_cell", "raw_background")
    groups = [
        (f"D{index:02d}", condition, label)
        for index in range(1, 4)
        for condition in ("Healthy", "Disease")
        for label in ("TypeA", "TypeB")
    ]
    assignments = [groups[index % len(groups)] for index in range(int(selected_mask.sum()))]
    columns = ("synthetic_donor", "synthetic_condition", "synthetic_label")
    for column, values in zip(columns, zip(*assignments, strict=True), strict=True):
        obs[column] = "background"
        obs.loc[selected_mask, column] = list(values)
    obs["original_total_counts"] = _row_sums(cast(Any, adata.X))
    obs["original_n_genes_by_counts"] = _row_nnz(cast(Any, adata.X))


def _assert_pair_contract(filtered: ad.AnnData, raw: ad.AnnData) -> None:
    if not filtered.var_names.equals(raw.var_names):
        raise ContractError("Filtered/raw pair has mismatched feature axes")
    filtered_obs = cast(pd.DataFrame, filtered.obs)
    raw_obs = cast(pd.DataFrame, raw.obs)
    raw_selected = raw_obs.loc[raw_obs["is_selected_filtered_cell"], "source_barcode"]
    if set(filtered_obs["source_barcode"]) != set(raw_selected):
        raise ContractError("Filtered/raw pair has mismatched selected barcodes")
    raw_background = raw_obs.loc[~raw_obs["is_selected_filtered_cell"], "source_role"]
    if not (raw_background == "raw_background").all():
        raise ContractError("Raw background droplets lost their explicit source role")
    aliases = ("synthetic_donor", "synthetic_condition", "synthetic_label")
    raw_selected_design = raw_obs.loc[raw_obs["is_selected_filtered_cell"]].set_index(
        "source_barcode"
    )[list(aliases)]
    filtered_design = filtered_obs.set_index("source_barcode")[list(aliases)]
    if not raw_selected_design.equals(filtered_design.loc[raw_selected_design.index]):
        raise ContractError("Raw selected cells lost their synthetic design metadata")


def _assert_individual_source_metadata(adata: ad.AnnData, sample: str) -> None:
    obs = cast(pd.DataFrame, adata.obs)
    reserved = {"batch", "donor", "condition", "label"}
    present = reserved & set(obs.columns)
    if present:
        raise ContractError(
            f"Individual source {sample} has reserved columns: {', '.join(sorted(present))}"
        )
    required = {
        "sample",
        "source_sample",
        "source_barcode",
        "synthetic_donor",
        "synthetic_condition",
        "synthetic_label",
        "source_role",
    }
    missing = required - set(obs.columns)
    if missing:
        raise ContractError(
            f"Individual source {sample} is missing metadata: {', '.join(sorted(missing))}"
        )
    expected_prefix = f"{sample}::"
    malformed = any(
        not name.startswith(expected_prefix) or name.count("::") != 1 for name in adata.obs_names
    )
    if malformed:
        raise ContractError(f"Individual source {sample} has malformed prefixed barcodes")


def _materialise_combined_design(adata: ad.AnnData) -> None:
    obs = cast(pd.DataFrame, adata.obs)
    obs["batch"] = obs["source_sample"]
    obs["donor"] = obs["synthetic_donor"]
    obs["condition"] = obs["synthetic_condition"]
    obs["label"] = obs["synthetic_label"]


def _write_csv(adata: ad.AnnData, path: Path) -> None:
    matrix: Any = adata[:CSV_CELLS, :CSV_GENES].X
    values = matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)
    frame = pd.DataFrame(
        values,
        index=adata.obs_names[:CSV_CELLS],
        columns=adata.var_names[:CSV_GENES],
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, float_format="%.0f")


def _write_10x_h5(adata: ad.AnnData, path: Path) -> None:
    x: Any = adata.X
    matrix = x.tocsc() if sparse.issparse(x) else sparse.csc_matrix(x)
    tenx_matrix = matrix.transpose().tocsc()
    if not np.allclose(matrix.data, np.rint(matrix.data)):
        raise ContractError("10x representation requires integer count values")
    path.parent.mkdir(parents=True, exist_ok=True)
    strings = h5py.string_dtype(encoding="utf-8")
    with h5py.File(path, "w") as handle:
        group = handle.create_group("matrix")
        group.create_dataset(
            "data", data=np.rint(tenx_matrix.data).astype(np.int32), compression="gzip"
        )
        group.create_dataset(
            "indices", data=tenx_matrix.indices.astype(np.int32), compression="gzip"
        )
        group.create_dataset(
            "indptr", data=tenx_matrix.indptr.astype(np.int32), compression="gzip"
        )
        group.create_dataset("shape", data=np.asarray(tenx_matrix.shape, dtype=np.int64))
        group.create_dataset(
            "barcodes", data=adata.obs_names.astype(str).to_numpy(dtype=object), dtype=strings
        )
        features = group.create_group("features")
        gene_names = adata.var_names.astype(str).to_numpy(dtype=object)
        features.create_dataset("id", data=gene_names, dtype=strings)
        features.create_dataset("name", data=gene_names, dtype=strings)
        features.create_dataset(
            "feature_type",
            data=np.repeat("Gene Expression", adata.n_vars).astype(object),
            dtype=strings,
        )
        features.create_dataset(
            "genome", data=np.repeat("unknown", adata.n_vars).astype(object), dtype=strings
        )


def _build_rds_sources(context: BuildContext, sources: dict[str, Path]) -> list[dict[str, Any]]:
    runtime = context.contracts["sources"]["pbmc"]["r_runtime"]
    helper = Path(__file__).with_name("source_formats.R")
    gene_panel = context.out_dir / "pbmc" / "combined_filtered_matrix.h5ad"
    if not gene_panel.is_file():
        raise ContractError("RDS conversion requires the selected combined gene panel")
    observations: list[dict[str, Any]] = []
    for sample, format_name in RDS_SAMPLES:
        suffix = "sce.rds" if format_name == "sce" else "seurat.rds"
        filtered_name = f"{sample}_filtered_matrix.{suffix}"
        raw_name = f"{sample}_raw_matrix.{suffix}"
        filtered_output = context.out_dir / "pbmc" / filtered_name
        raw_output = context.out_dir / "pbmc" / raw_name
        _run_rds_conversion(
            helper,
            runtime,
            sources[filtered_name],
            sources[raw_name],
            filtered_output,
            raw_output,
            sample,
            format_name,
            context.cache_dir / "source-work" / "rds-conversions",
            gene_panel,
        )
        observations.append(
            {
                "sample": sample,
                "format": format_name,
                "filtered_bytes": filtered_output.stat().st_size,
                "raw_bytes": raw_output.stat().st_size,
                "container": runtime["container"],
            }
        )
    return observations


def _run_rds_conversion(
    helper: Path,
    runtime: dict[str, str],
    filtered_input: Path,
    raw_input: Path,
    filtered_output: Path,
    raw_output: Path,
    sample: str,
    format_name: str,
    work_root: Path,
    gene_panel: Path,
) -> None:
    filtered_output.parent.mkdir(parents=True, exist_ok=True)
    # Apptainer can bind the shared cache, keeping producer scratch off host /tmp.
    work_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="conversion-", dir=work_root) as temporary:
        work = Path(temporary)
        for source in (helper, filtered_input, raw_input):
            shutil.copy2(source, work / source.name)
        panel_adata = ad.read_h5ad(gene_panel, backed="r")
        try:
            (work / "selected-genes.txt").write_text(
                "\n".join(panel_adata.var_names.astype(str)) + "\n", encoding="utf-8"
            )
        finally:
            panel_adata.file.close()
        command = [
            "apptainer",
            "exec",
            "--containall",
            "--home",
            f"{work}:/work",
            "--bind",
            f"{work}:/work",
            "--pwd",
            "/work",
            _apptainer_image(runtime["container"]),
            "Rscript",
            "/work/source_formats.R",
            format_name,
            f"/work/{filtered_input.name}",
            f"/work/{raw_input.name}",
            "/work/filtered.rds",
            "/work/raw.rds",
            sample,
            "384",
            "3000",
            str(RAW_BACKGROUND_CELLS),
            "/work/selected-genes.txt",
        ]
        try:
            subprocess.run(command, check=True, text=True, capture_output=True)
        except FileNotFoundError as error:
            raise ContractError("Apptainer is required for pinned RDS source conversion") from error
        except subprocess.CalledProcessError as error:
            raise ContractError(
                f"Pinned RDS conversion failed for {sample}: {error.stderr.strip()}"
            ) from error
        shutil.copy2(work / "filtered.rds", filtered_output)
        shutil.copy2(work / "raw.rds", raw_output)


def _enforce_size_limit(context: BuildContext, produced: list[str]) -> None:
    hard_limit = context.contracts["size_policy"]["hard_limit_bytes"]
    too_large = [
        f"{path} ({(context.out_dir / path).stat().st_size} bytes)"
        for path in produced
        if (context.out_dir / path).stat().st_size >= hard_limit
    ]
    if too_large:
        raise ContractError(f"Source stage exceeded the hard size limit: {', '.join(too_large)}")


def _verify_production_readers(context: BuildContext) -> list[dict[str, Any]]:
    readers = context.contracts["sources"]["pbmc"]["reader_runtimes"]
    templates = {
        "csv": context.pipeline_dir / "modules/local/adata/readcsv/templates/readcsv.py",
        "tenx": context.pipeline_dir / "modules/local/scanpy/readh5/templates/readh5.py",
        "rds": context.pipeline_dir / "modules/local/adata/readrds/templates/readrds.py",
    }
    work_root = context.cache_dir / "source-work" / "reader-probes"
    work_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="probe-", dir=work_root) as temporary:
        work = Path(temporary)
        outputs: list[dict[str, Any]] = []
        csv = context.out_dir / "pbmc/SRR28679756_filtered_matrix.csv"
        _copy_reader_template(
            templates["csv"],
            work / "reader_readcsv.py",
            {"csv": "/work/input.csv", "prefix": "/work/csv"},
        )
        shutil.copy2(csv, work / "input.csv")
        _run_container_python(work, readers["csv"]["container"], "/work/reader_readcsv.py")
        outputs.append(
            _reader_output("csv", readers["csv"]["container"], work / "csv.h5ad", (64, 128))
        )

        tenx = context.out_dir / "misc/filtered_feature_bc_matrix.h5"
        _copy_reader_template(
            templates["tenx"], work / "tenx.py", {"h5": "/work/input.h5", "prefix": "/work/tenx"}
        )
        shutil.copy2(tenx, work / "input.h5")
        _run_container_python(work, readers["tenx"]["container"], "/work/tenx.py")
        outputs.append(
            _reader_output("tenx", readers["tenx"]["container"], work / "tenx.h5ad", (384, 3000))
        )

        for sample, format_name in RDS_SAMPLES:
            suffix = "sce.rds" if format_name == "sce" else "seurat.rds"
            source = context.out_dir / "pbmc" / f"{sample}_filtered_matrix.{suffix}"
            name = f"{sample}_{format_name}"
            _copy_reader_template(
                templates["rds"],
                work / f"{name}.py",
                {"rds": f"/work/{name}.rds", "prefix": f"/work/{name}"},
            )
            shutil.copy2(source, work / f"{name}.rds")
            _run_container_python(work, readers["rds"]["container"], f"/work/{name}.py")
            outputs.append(
                _reader_output(
                    name,
                    readers["rds"]["container"],
                    work / f"{name}.h5ad",
                    (384, 3000),
                    forbidden_obs={"batch", "donor", "condition", "label"},
                )
            )
    return outputs


def _copy_reader_template(template: Path, destination: Path, values: dict[str, str]) -> None:
    script = template.read_text(encoding="utf-8")
    replacements = {**values, "task.process": "SOURCE_READER_PROBE"}
    for key, value in replacements.items():
        script = script.replace(f"${{{key}}}", value)
    destination.write_text(script, encoding="utf-8")


def _run_container_python(work: Path, image: str, script: str) -> None:
    command = [
        "apptainer",
        "exec",
        "--containall",
        "--home",
        f"{work}:/work",
        "--env",
        "PATH=/opt/wave/.pixi/envs/default/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin",
        "--env",
        "MPLCONFIGDIR=/work/.mplconfig,NUMBA_CACHE_DIR=/work/.numba",
        "--bind",
        f"{work}:/work",
        "--pwd",
        "/work",
        _apptainer_image(image),
        "python",
        script,
    ]
    try:
        subprocess.run(command, check=True, text=True, capture_output=True)
    except FileNotFoundError as error:
        raise ContractError("Apptainer is required for production reader probes") from error
    except subprocess.CalledProcessError as error:
        raise ContractError(f"Production reader probe failed: {error.stderr.strip()}") from error


def _reader_output(
    name: str,
    container: str,
    path: Path,
    shape: tuple[int, int],
    forbidden_obs: set[str] | None = None,
    required_var: set[str] | None = None,
) -> dict[str, Any]:
    output = ad.read_h5ad(path)
    if output.shape != shape:
        raise ContractError(f"Production reader '{name}' returned {output.shape}, expected {shape}")
    present = (forbidden_obs or set()) & set(cast(pd.DataFrame, output.obs).columns)
    if present:
        raise ContractError(
            f"Production reader '{name}' restored reserved columns: {', '.join(sorted(present))}"
        )
    missing_var = (required_var or set()) - set(cast(pd.DataFrame, output.var).columns)
    if missing_var:
        raise ContractError(
            f"Production reader '{name}' omitted gene symbols: {', '.join(sorted(missing_var))}"
        )
    return {"reader": name, "container": container, "shape": list(output.shape)}


def _apptainer_image(image: str) -> str:
    """Return an explicit Docker transport URI for a pinned registry image."""
    return image if "://" in image else f"docker://{image}"


def _h5ad_observations(out_dir: Path, produced: list[str]) -> dict[str, Any]:
    observations: dict[str, Any] = {}
    for relative in produced:
        if not relative.endswith(".h5ad"):
            continue
        path = out_dir / relative
        violations = h5ad_compression_violations(path)
        if violations:
            raise ContractError(
                f"Source H5AD has uncompressed arrays: {relative}: {', '.join(violations)}"
            )
        adata = ad.read_h5ad(path)
        observations[relative] = {
            "shape": list(adata.shape),
            "x_format": "sparse" if sparse.issparse(adata.X) else "dense",
            "x_dtype": str(cast(Any, adata.X).dtype),
            "bytes": path.stat().st_size,
        }
    return observations


def _row_sums(matrix: Any) -> np.ndarray:
    return np.asarray(matrix.sum(axis=1)).ravel().astype(np.int64)


def _row_nnz(matrix: Any) -> np.ndarray:
    return np.asarray((matrix > 0).sum(axis=1)).ravel().astype(np.int64)


def _variance(matrix: Any) -> np.ndarray:
    matrix = matrix.tocsr().astype(np.float64)
    mean = np.asarray(matrix.mean(axis=0)).ravel()
    mean_squared = np.asarray(matrix.multiply(matrix).mean(axis=0)).ravel()
    return mean_squared - mean**2


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
