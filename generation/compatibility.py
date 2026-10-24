"""Materialise a separately checksummed compatibility candidate without changing counts."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import anndata as ad
import numpy as np

from .core import ContractError, sha256_file, write_h5ad


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = args.input.resolve()
    destination = args.output.resolve()
    if destination.exists() or destination.is_relative_to(source):
        raise ContractError("Use a new output directory outside the input candidate")
    destination.mkdir(parents=True)
    records = []
    for path in sorted(source.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(source)
        output = destination / relative
        if relative.parts[0] == "reports":
            output = destination / "reports" / "upstream" / relative.name
        output.parent.mkdir(parents=True, exist_ok=True)
        if path.suffix == ".h5ad":
            original = ad.read_h5ad(path)
            repaired = original.copy()
            if relative.parts[0] == "pbmc":
                repaired.var["symbols"] = repaired.var_names.to_numpy(dtype=object)
            write_h5ad(repaired, output)
            restored = ad.read_h5ad(output)
            if original.obs_names.tolist() != restored.obs_names.tolist():
                raise ContractError(f"Cell identities changed: {relative}")
            if original.var_names.tolist() != restored.var_names.tolist():
                raise ContractError(f"Feature identities changed: {relative}")
            for name in ["X", *original.layers]:
                before = original.X if name == "X" else original.layers[name]
                after = restored.X if name == "X" else restored.layers[name]
                difference = before != after
                unequal = difference.nnz if hasattr(difference, "nnz") else np.any(difference)
                if unequal:
                    raise ContractError(f"Expression changed in {relative}: {name}")
        else:
            shutil.copy2(path, output)
        records.append({
            "path": relative.as_posix(), "input_sha256": sha256_file(path),
            "output_sha256": sha256_file(output), "bytes": output.stat().st_size,
        })
    report = destination / "reports" / "compatibility.json"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(json.dumps({
        "input": str(source),
        "transformation": (
            "Legacy string encoding; explicit PBMC symbols; unchanged expression and axes"
        ),
        "artefacts": records,
        "publication_status": "Acceptance candidate; complete stage provenance remains required",
    }, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"candidate": str(destination), "files": len(records), "status": "passed"}))


if __name__ == "__main__":
    main()
