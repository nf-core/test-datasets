"""Generate the matching reference triplet with the production pipeline."""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from .core import BuildContext, ContractError, StageResult, sha256_file, write_h5ad

CONFIG = """params {
    ambient_correction      = 'none'
    doublet_detection       = 'none'
    cell_cycle_scoring      = false
    celltypist_model        = null
    celldex_reference       = null
    cytetype_study_context  = ''
    scib                    = false
    skip_liana              = true
    skip_qc_report          = true
    integration_methods     = 'scvi,symphony'
    normalization_method    = 'log1p'
    feature_selection       = 'hvgs'
    integration_n_features  = 500
    scvi_max_epochs         = 2
    scvi_n_hidden           = 16
    scvi_n_layers           = 1
    scvi_n_latent           = 5
    clustering_resolutions  = '0.5'
    cluster_global          = false
    cluster_per_label       = false
    analysis_plan           = null
    de_methods              = ''
    sample_n                = 128
}
tower.enabled = false
process.beforeScript = 'source /etc/profile.d/modules.sh; module load apptainer'
"""

INPUTS = (
    ("SRR28679756", "SRR28679756_filtered_matrix.h5ad"),
    ("SRR28679757", "SRR28679757_filtered_matrix.sce.rds"),
    ("SRR28679758", "SRR28679758_filtered_matrix.seurat.rds"),
)


def build(context: BuildContext) -> StageResult:
    """Train all reference artefacts together from rebuilt source representations."""
    import anndata as ad

    sources = [context.out_dir / "pbmc" / name for _, name in INPUTS]
    missing = [str(path) for path in sources if not path.is_file()]
    if missing:
        raise ContractError(f"Reference generation requires source-stage outputs: {missing}")
    nextflow = shutil.which("nextflow")
    if nextflow is None:
        raise ContractError("Provision Nextflow on PATH before reference generation")
    environment = os.environ.copy()
    environment["NXF_VER"] = context.contracts["pipeline"]["nextflow"]
    environment["NXF_ANSI_LOG"] = "false"
    version = subprocess.check_output(
        [nextflow, "-version"], env=environment, text=True, stderr=subprocess.STDOUT
    )
    if f"version {environment['NXF_VER']} " not in version:
        raise ContractError(f"Nextflow does not match the locked runtime: {version.strip()}")
    # Slurm tasks must see the same shared paths as the Nextflow driver.
    scratch_root = Path(os.environ.get(
        "SCDOWNSTREAM_CONTAINER_WORK_ROOT", str(context.cache_dir / "reference-work")
    ))
    scratch_root.mkdir(parents=True, exist_ok=True)
    logs = context.cache_dir / "reference-runs" / context.started_at.replace(":", "-")
    logs.mkdir(parents=True, exist_ok=False)
    source_hashes = {
        path.relative_to(context.out_dir).as_posix(): sha256_file(path) for path in sources
    }
    with tempfile.TemporaryDirectory(
        prefix="scdownstream-reference-", dir=scratch_root
    ) as temporary:
        root = Path(temporary)
        pipeline = root / "pipeline"
        shutil.copytree(
            context.pipeline_dir,
            pipeline,
            ignore=shutil.ignore_patterns(".git", ".nf-test", ".nextflow", "work", "results"),
        )
        inputs = root / "inputs"
        inputs.mkdir()
        sheet = root / "samplesheet.csv"
        with sheet.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["sample", "filtered", "batch_col"])
            for (sample, _), source in zip(INPUTS, sources, strict=True):
                destination = inputs / source.name
                shutil.copy2(source, destination)
                writer.writerow([sample, str(destination), "sample"])
        config = root / "generate.config"
        config.write_text(CONFIG, encoding="utf-8")
        results = root / "results"
        command = [
            nextflow,
            "run",
            str(pipeline / "main.nf"),
            "-profile",
            "test,apptainer,daisybio,keep_work",
            "-c",
            str(config),
            "-work-dir",
            str(root / "work"),
            "--input",
            str(sheet),
            "--outdir",
            str(results),
        ]
        (logs / "generate.config").write_text(CONFIG, encoding="utf-8")
        (logs / "command.txt").write_text("\n".join(command) + "\n", encoding="utf-8")
        try:
            with (logs / "nextflow.stdout.log").open("w", encoding="utf-8") as output:
                subprocess.run(
                    command,
                    cwd=root,
                    env=environment,
                    stdout=output,
                    stderr=subprocess.STDOUT,
                    check=True,
                )
        finally:
            if (root / ".nextflow.log").is_file():
                shutil.copy2(root / ".nextflow.log", logs / "nextflow.log")
            if (results / "pipeline_info").is_dir():
                shutil.copytree(results / "pipeline_info", logs / "pipeline_info")
        merged = ad.read_h5ad(results / "09_finalized" / "merged.h5ad")
        if "raw_counts" not in merged.layers:
            raise ContractError("Production reference did not retain raw_counts")
        if set(merged.obs["sample"].astype(str)) != {sample for sample, _ in INPUTS}:
            raise ContractError("Production reference lost a source sample")
        reference = ad.read_h5ad(results / "03_combine/integrate/symphony/symphony_reference.h5ad")
        out = context.out_dir / "extension_base"
        write_h5ad(merged, out / "merged.h5ad")
        write_h5ad(reference, out / "symphony_reference.h5ad")
        shutil.copy2(results / "03_combine/integrate/scvi/scvi_model/model.pt", out / "model.pt")
        observations = {
            "n_obs": merged.n_obs,
            "n_vars": merged.n_vars,
            "symphony_n_vars": reference.n_vars,
            "mapping_and_extension": "Require separate production nf-test acceptance runs",
        }
    produced = [
        f"extension_base/{name}" for name in ("merged.h5ad", "model.pt", "symphony_reference.h5ad")
    ]
    return StageResult(
        stage="reference",
        produced=produced,
        provenance={
            "pipeline_revision": context.pipeline_revision,
            "nextflow": environment["NXF_VER"],
            "source_sha256": source_hashes,
            "model_architecture": {"hidden": 16, "layers": 1, "latent": 5, "epochs": 2},
            "source_selection_seed": context.seed,
            "scvi_seed": 0,
            "artefact_sha256": {name: sha256_file(context.out_dir / name) for name in produced},
        },
        dependency_fingerprint=context.dependency_fingerprint,
        observations=observations,
    )
