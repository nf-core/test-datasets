#!/usr/bin/env python3
"""Derive canonical PCA fixtures from canonical source-stage count H5ADs."""

from __future__ import annotations

import argparse
from pathlib import Path

from generation.analysis import PCA_INPUTS, derive_pca_files


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-root", type=Path, required=True, help="candidate root containing pbmc/"
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True, help="external directory for PCA H5ADs"
    )
    args = parser.parse_args()
    source_root = args.source_root.resolve()
    output_dir = args.output_dir.resolve()
    if output_dir == source_root / "pbmc":
        parser.error("Use an output directory separate from source-stage count fixtures")
    inputs = {output: source_root / "pbmc" / source for output, source in PCA_INPUTS.items()}
    missing = [str(path) for path in inputs.values() if not path.is_file()]
    if missing:
        parser.error("Missing canonical source-stage inputs: " + ", ".join(missing))
    outputs = {name: output_dir / name for name in PCA_INPUTS}
    work_root = output_dir.parent / "pca-fixture-work"
    observations = derive_pca_files(inputs, outputs, work_root)
    for name, observation in observations.items():
        print(f"{name}: {observation['shape']} ({observation['bytes']} bytes, gzip)")


if __name__ == "__main__":
    main()
