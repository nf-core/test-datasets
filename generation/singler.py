"""Build compact authentic SingleR reference archives."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import re
import shutil
import subprocess
import tarfile
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.error import URLError
from urllib.parse import urlparse
from urllib.request import urlopen

import anndata as ad
import h5py

from .core import BuildContext, ContractError, StageResult, sha256_file

R_CONTAINER = (
    "community.wave.seqera.io/library/bioconductor-anndatar_bioconductor-hdf5array_"
    "bioconductor-rhdf5_bioconductor-singlecellexperiment_pruned@sha256:"
    "0bc8f1a815994cbc7527a4d25a7cbda31efd071879a3f14c3a8e712e6b87dd8f"
)
DATASET_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE_REVISION = "e3a7f43eb802a090affac918026d2ba5dce8fcd5"
MINIMUM_SHARED_GENES = 1_000
TARGET_GENES = 1_000

SOURCE_ARCHIVES: tuple[dict[str, Any], ...] = (
    {
        "name": "celldex_hpca__2024-02-26_h5_se.tar",
        "uri": (
            "https://raw.githubusercontent.com/nf-core/test-datasets/"
            f"{DEFAULT_SOURCE_REVISION}/singleR/celldex_hpca__2024-02-26_h5_se.tar"
        ),
        "sha256": "85d2fa25cb81fcae72594b64e55b4ceccf52f1639100d22d3c0d428d63b7c007",
        "reference": "hpca",
        "label_field": "label.fine",
        "stratify_field": "label.fine",
        "profile_count": 128,
        "minimum_label_categories": 128,
    },
    {
        "name": "celldex_monaco_immune__2024-02-26_h5_se.tar",
        "uri": (
            "https://raw.githubusercontent.com/nf-core/test-datasets/"
            f"{DEFAULT_SOURCE_REVISION}/singleR/celldex_monaco_immune__2024-02-26_h5_se.tar"
        ),
        "sha256": "3d3e207cd6741c4e68cd8165d273eee1bdad183754fad44c62423af132d644e8",
        "reference": "monaco_immune",
        "label_field": "label.main",
        "stratify_field": "label.fine",
        "profile_count": 114,
        "minimum_label_categories": 10,
    },
)


@dataclass(frozen=True)
class ArchiveSpec:
    """Pinned source and retained-profile contract for one reference."""

    name: str
    uri: str
    sha256: str
    reference: str
    label_field: str
    stratify_field: str
    profile_count: int
    minimum_label_categories: int

    @property
    def stem(self) -> str:
        return self.name.removesuffix(".tar")


def build(context: BuildContext) -> StageResult:
    """Subset pinned celldex archives against the produced PBMC gene panel."""
    single_r = context.out_dir / "singleR"
    source_specs = _source_specs(context)
    source_revision = _source_revision(source_specs)
    singler_contract = context.contracts.get("sources", {}).get("singler", {})
    target_genes = int(singler_contract.get("target_genes", TARGET_GENES))
    minimum_shared_genes = int(
        singler_contract.get("minimum_shared_genes", MINIMUM_SHARED_GENES)
    )
    if target_genes <= 0 or minimum_shared_genes < target_genes:
        raise ContractError("SingleR gene counts must be positive and overlap must cover selection")
    panel_path = _gene_panel_path(context)
    panel_sha256 = sha256_file(panel_path)
    panel_genes = _read_gene_panel(panel_path, minimum_shared_genes)
    references_sheet = DATASET_ROOT / "singleR" / "references.csv"
    _read_references_sheet(references_sheet)
    target_bytes = context.contracts["size_policy"]["target_bytes"]
    hard_limit = context.contracts["size_policy"]["hard_limit_bytes"]

    cache_source_dir = context.cache_dir / "sources" / "singler"
    cache_source_dir.mkdir(parents=True, exist_ok=True)
    image_id = _apptainer_image(R_CONTAINER)
    producer = Path(__file__).with_name("singler_formats.R")
    if not producer.is_file():
        raise ContractError(f"SingleR R producer is missing: {producer}")

    work_root = context.cache_dir / "singler-work"
    work_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="build-", dir=work_root) as temporary:
        scratch = Path(temporary)
        producer_root = scratch / "producer"
        producer_root.mkdir()
        shutil.copy2(producer, producer_root / producer.name)
        panel_path_local = producer_root / "query-genes.txt"
        panel_path_local.write_text("\n".join(panel_genes) + "\n", encoding="utf-8")
        query_path_local = producer_root / panel_path.name
        shutil.copyfile(panel_path, query_path_local)
        panel_hash = _sha256_bytes(panel_path_local.read_bytes())

        summaries: list[dict[str, Any]] = []
        archive_records: list[dict[str, Any]] = []
        for spec in source_specs:
            source_archive = _resolve_source_archive(context, spec, cache_source_dir)
            local_archive = producer_root / spec.name
            shutil.copyfile(source_archive, local_archive)
            output_dir = producer_root / f"output-{spec.reference}"
            output_dir.mkdir()
            result_path = producer_root / f"{spec.reference}-probe.json"
            _run_r_producer(
                local_archive,
                panel_path_local,
                query_path_local,
                output_dir,
                result_path,
                spec,
                panel_hash,
                target_genes,
            )
            summary = json.loads(result_path.read_text(encoding="utf-8"))
            if summary["shared_gene_count"] < minimum_shared_genes:
                raise ContractError(
                    f"{spec.name} has only {summary['shared_gene_count']} shared genes"
                )
            if summary["gene_count"] != target_genes:
                raise ContractError(f"{spec.name} did not retain {target_genes} genes")
            if summary["label_categories"] < spec.minimum_label_categories:
                raise ContractError(
                    f"{spec.name} retains {summary['label_categories']} "
                    f"{spec.label_field} categories, expected at least "
                    f"{spec.minimum_label_categories}"
                )

            staged_ref = output_dir / spec.stem
            archive_path = single_r / spec.name
            archive_path.parent.mkdir(parents=True, exist_ok=True)
            _create_reproducible_tar(staged_ref, archive_path)
            compressed_datasets = _validate_archive(archive_path, spec)
            archive_size = archive_path.stat().st_size
            if archive_size >= hard_limit:
                raise ContractError(
                    f"SingleR archive exceeds the hard limit: {spec.name} ({archive_size} bytes)"
                )
            archive_record = {
                "path": archive_path.relative_to(context.out_dir).as_posix(),
                "bytes": archive_size,
                "sha256": sha256_file(archive_path),
                "source_path": str(source_archive),
                "source_sha256": spec.sha256,
                "source_uri": spec.uri,
                "gzip_assay_datasets": compressed_datasets,
                **summary,
            }
            archive_records.append(archive_record)
            summaries.append(summary)

        generated_sheet = single_r / "references.csv"
        generated_sheet.write_bytes(references_sheet.read_bytes())
        produced = [record["path"] for record in archive_records]
        produced.append(generated_sheet.relative_to(context.out_dir).as_posix())
        target_reviews = [
            {"path": record["path"], "bytes": record["bytes"]}
            for record in archive_records
            if record["bytes"] >= target_bytes
        ]
        if target_reviews:
            raise ContractError(
                "SingleR archives exceed the target size: "
                + ", ".join(f"{item['path']}={item['bytes']}" for item in target_reviews)
            )

    return StageResult(
        stage="singler",
        produced=produced,
        provenance={
            "source_archives": archive_records,
            "source_revision": source_revision,
            "gene_panel": {
                "path": panel_path.relative_to(context.out_dir).as_posix(),
                "file_sha256": panel_sha256,
                "ordered_gene_list_sha256": panel_hash,
                "source_genes": len(panel_genes),
                "selected_genes_per_reference": target_genes,
                "selection": "first shared genes in source var_names order",
            },
            "references_sheet_sha256": sha256_file(references_sheet),
            "container": R_CONTAINER,
            "container_image_id": image_id,
            "runtime": summaries[0]["runtime"],
        },
        dependency_fingerprint=context.dependency_fingerprint,
        observations={
            "references": archive_records,
            "query_probes": [summary["inference"] for summary in summaries],
            "target_size_reviews": target_reviews,
        },
    )


def _source_specs(context: BuildContext) -> list[ArchiveSpec]:
    source_contract = context.contracts.get("sources", {}).get("singler", {})
    configured = source_contract.get("archives", source_contract.get("uris", []))
    entries = configured or SOURCE_ARCHIVES
    specs: list[ArchiveSpec] = []
    for entry in entries:
        if not isinstance(entry, dict):
            raise ContractError("SingleR source entries must be objects with URI and SHA256")
        name = entry.get("name") or Path(entry.get("path", "")).name
        fallback = next(
            (candidate for candidate in SOURCE_ARCHIVES if candidate["name"] == name), None
        )
        if fallback is None:
            raise ContractError(f"Unsupported SingleR archive name in contract: {name!r}")
        spec = ArchiveSpec(
            name=name,
            uri=entry.get("uri", fallback["uri"]),
            sha256=entry.get("sha256", fallback["sha256"]),
            reference=entry.get("reference", fallback["reference"]),
            label_field=entry.get("label_field", fallback["label_field"]),
            stratify_field=entry.get("stratify_field", fallback["stratify_field"]),
            profile_count=int(entry.get("profile_count", fallback["profile_count"])),
            minimum_label_categories=int(
                entry.get("minimum_label_categories", fallback["minimum_label_categories"])
            ),
        )
        if len(spec.sha256) != 64 or any(
            character not in "0123456789abcdef" for character in spec.sha256
        ):
            raise ContractError(f"Invalid SHA256 for SingleR source {spec.name}")
        specs.append(spec)
    if {spec.reference for spec in specs} != {"hpca", "monaco_immune"}:
        raise ContractError("SingleR stage requires the pinned HPCA and Monaco references")
    if len({spec.name for spec in specs}) != len(specs):
        raise ContractError("SingleR source archive names must be unique")
    return specs


def _source_revision(specs: list[ArchiveSpec]) -> str | None:
    revisions: set[str] = set()
    for spec in specs:
        path_parts = urlparse(spec.uri).path.strip("/").split("/")
        for index, part in enumerate(path_parts[:-1]):
            if part == "test-datasets" and re.fullmatch(r"[0-9a-f]{40}", path_parts[index + 1]):
                revisions.add(path_parts[index + 1])
    if len(revisions) > 1:
        raise ContractError("SingleR source archives must use the same pinned dataset revision")
    return next(iter(revisions), None)


def _gene_panel_path(context: BuildContext) -> Path:
    singler_contract = context.contracts.get("sources", {}).get("singler", {})
    pbmc_contract = context.contracts.get("sources", {}).get("pbmc", {})
    relative = singler_contract.get("gene_panel_path") or pbmc_contract.get("gene_panel_path")
    if not relative:
        raise ContractError(
            "SingleR requires sources.singler.gene_panel_path or "
            "sources.pbmc.gene_panel_path to declare the produced query H5AD"
        )
    path = (context.out_dir / relative).resolve()
    if not path.is_relative_to(context.out_dir) or not path.is_file():
        raise ContractError(f"Declared SingleR gene-panel H5AD is missing: {path}")
    return path


def _read_gene_panel(path: Path, minimum_shared_genes: int) -> list[str]:
    try:
        matrix = ad.read_h5ad(path, backed="r")
        genes = [str(gene) for gene in matrix.var_names]
        matrix.file.close()
    except (OSError, ValueError) as error:
        raise ContractError(f"Cannot read the SingleR query gene panel {path}: {error}") from error
    if not genes or any(not gene or gene == "nan" for gene in genes):
        raise ContractError(f"SingleR query gene panel contains empty names: {path}")
    if len(set(genes)) != len(genes):
        raise ContractError(f"SingleR query gene panel contains duplicate symbols: {path}")
    if len(genes) < minimum_shared_genes:
        raise ContractError(
            f"SingleR query gene panel has {len(genes)} genes, below {minimum_shared_genes}"
        )
    return genes


def _resolve_source_archive(
    context: BuildContext, spec: ArchiveSpec, cache_source_dir: Path
) -> Path:
    cached = cache_source_dir / spec.name
    local = DATASET_ROOT / "singleR" / spec.name
    for candidate in (cached, local):
        if candidate.is_file() and sha256_file(candidate) == spec.sha256:
            return candidate
    if not spec.uri:
        raise ContractError(f"No verified source archive or URI is available for {spec.name}")
    cached.parent.mkdir(parents=True, exist_ok=True)
    temporary = cached.with_suffix(cached.suffix + ".download")
    digest = hashlib.sha256()
    try:
        with urlopen(spec.uri, timeout=120) as response, temporary.open("wb") as output:
            while block := response.read(1024 * 1024):
                digest.update(block)
                output.write(block)
    except (OSError, URLError) as error:
        temporary.unlink(missing_ok=True)
        raise ContractError(f"Cannot download pinned SingleR source {spec.uri}: {error}") from error
    if digest.hexdigest() != spec.sha256:
        temporary.unlink(missing_ok=True)
        raise ContractError(
            f"Downloaded SingleR source checksum mismatch for {spec.name}: {digest.hexdigest()}"
        )
    temporary.replace(cached)
    return cached


def _read_references_sheet(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise ContractError(f"SingleR reference projection is missing: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or set(rows[0]) != {"id", "label", "reference", "version"}:
        raise ContractError(f"SingleR reference projection has an unexpected schema: {path}")
    ids = {row["id"] for row in rows}
    required_ids = {"hpca_celldex", "immune_celldex", "hpca_direct", "immune_direct"}
    if not required_ids.issubset(ids):
        raise ContractError(f"SingleR reference projection is missing rows: {path}")
    for row in rows:
        if row["id"] == "hpca_celldex" and (row["label"], row["reference"]) != (
            "label.main",
            "hpca",
        ):
            raise ContractError(f"Unexpected HPCA source-sheet projection in {path}")
        if row["id"] == "immune_celldex" and (
            row["label"],
            row["reference"],
        ) != ("label.fine", "monaco_immune"):
            raise ContractError(f"Unexpected Monaco source-sheet projection in {path}")
    return rows


def _validate_archive(archive_path: Path, spec: ArchiveSpec) -> int:
    expected_stem = spec.stem
    with tarfile.open(archive_path) as archive:
        members = set(archive.getnames())
        required_members = {
            expected_stem,
            f"{expected_stem}/assays.h5",
            f"{expected_stem}/se.rds",
        }
        if not required_members.issubset(members):
            raise ContractError(f"SingleR archive has an invalid stem layout: {archive_path}")
        assay_file = archive.extractfile(f"{expected_stem}/assays.h5")
        if assay_file is None:
            raise ContractError(f"SingleR archive has no assay file: {archive_path}")
        with h5py.File(io.BytesIO(assay_file.read()), "r") as handle:
            datasets: list[h5py.Dataset] = []

            def inspect(_name: str, node: h5py.Dataset | h5py.Group) -> None:
                if isinstance(node, h5py.Dataset):
                    datasets.append(node)

            handle.visititems(inspect)
            violations = [
                dataset.name
                for dataset in datasets
                if dataset.size > 0 and dataset.compression != "gzip"
            ]
            if not datasets or violations:
                raise ContractError(
                    f"SingleR assay datasets are missing gzip compression: "
                    f"{archive_path}: {violations}"
                )
            return len(datasets)


def _run_r_producer(
    archive: Path,
    panel_path: Path,
    query_h5ad: Path,
    output_dir: Path,
    result_path: Path,
    spec: ArchiveSpec,
    panel_hash: str,
    target_genes: int,
) -> None:
    # --containall provides a small isolated /tmp, unsuitable for source archives.
    (archive.parent / "tmp").mkdir(exist_ok=True)
    arguments = [
        "apptainer",
        "exec",
        "--containall",
        "--home",
        f"{archive.parent}:/work",
        "--bind",
        f"{archive.parent}:/work",
        "--env",
        "TMPDIR=/work/tmp",
        "--pwd",
        "/work",
        _apptainer_image(R_CONTAINER),
        "Rscript",
        "/work/singler_formats.R",
        f"/work/{archive.name}",
        f"/work/{panel_path.name}",
        f"/work/{query_h5ad.name}",
        f"/work/{output_dir.name}",
        f"/work/{result_path.name}",
        spec.stem,
        spec.label_field,
        spec.stratify_field,
        str(spec.profile_count),
        str(spec.minimum_label_categories),
        str(target_genes),
        panel_hash,
    ]
    try:
        subprocess.run(arguments, check=True, text=True, capture_output=True)
    except FileNotFoundError as error:
        raise ContractError("Apptainer is required to build SingleR references") from error
    except subprocess.CalledProcessError as error:
        message = error.stderr.strip()[-4000:]
        raise ContractError(f"SingleR R producer failed for {spec.name}: {message}") from error


def _apptainer_image(image: str) -> str:
    """Use an explicit Docker registry transport with the pinned image digest."""
    return image if "://" in image else f"docker://{image}"


def _create_reproducible_tar(source_dir: Path, destination: Path) -> None:
    if not source_dir.is_dir():
        raise ContractError(f"R producer did not create reference directory: {source_dir}")
    with tarfile.open(destination, mode="w", format=tarfile.USTAR_FORMAT) as archive:
        directory = tarfile.TarInfo(source_dir.name)
        directory.type = tarfile.DIRTYPE
        directory.mode = 0o755
        directory.mtime = 0
        archive.addfile(directory)
        for filename in ("assays.h5", "se.rds"):
            source_file = source_dir / filename
            if not source_file.is_file():
                raise ContractError(f"R producer omitted {source_file.name} from {source_dir}")
            info = archive.gettarinfo(str(source_file), arcname=f"{source_dir.name}/{filename}")
            info.uid = 0
            info.gid = 0
            info.uname = ""
            info.gname = ""
            info.mode = 0o644
            info.mtime = 0
            with source_file.open("rb") as handle:
                archive.addfile(info, handle)


def _sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()
