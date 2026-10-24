"""Generate small DE tables for the plotting regression scenarios."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

FILENAMES = [
    "two_groups_results.parquet",
    "two_groups_one_panel_results.parquet",
    "three_groups_results.parquet",
    "sparse_groups_results.parquet",
    "one_panel_results.parquet",
    "missing_pvalues_results.parquet",
    "raw_pvalues_results.parquet",
    "condition:leiden:11:wilcoxon_results.parquet",
]


def table(groups: list[str]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "gene": gene,
                "log2fc": fold_change,
                "pvalue": pvalue,
                "padj": pvalue * 2,
                "group": group,
                "contrast": f"{group} vs rest",
                "stratum": "condition=Healthy",
            }
            for group in groups
            for gene, fold_change, pvalue in [
                ("GENE1", 2.0, 0.001),
                ("GENE2", -1.5, 0.01),
                ("GENE3", 0.2, 0.5),
            ]
        ]
    )


def build_fixtures(output_dir: str | Path, seed: int = 42) -> list[str]:
    """Build eight deterministic Parquet edge cases without touching old outputs."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    two_groups = table(["B", "A"])
    two_groups_one_panel = two_groups.copy()
    two_groups_one_panel.loc[
        two_groups_one_panel["group"] == "A", ["pvalue", "padj"]
    ] = float("nan")
    three_groups = table(["C", "B", "A"])
    sparse_groups = three_groups.copy()
    sparse_groups.loc[sparse_groups["group"] == "C", ["pvalue", "padj"]] = float("nan")
    one_panel = sparse_groups.copy()
    one_panel.loc[one_panel["group"] == "A", ["pvalue", "padj"]] = float("nan")
    missing_pvalues = table(["A"])
    missing_pvalues[["pvalue", "padj"]] = float("nan")
    raw_pvalues = table(["A"])
    raw_pvalues["padj"] = float("nan")

    frames = {
        "two_groups": two_groups,
        "two_groups_one_panel": two_groups_one_panel,
        "three_groups": three_groups,
        "sparse_groups": sparse_groups,
        "one_panel": one_panel,
        "missing_pvalues": missing_pvalues,
        "raw_pvalues": raw_pvalues,
        "condition:leiden:11:wilcoxon": two_groups,
    }
    # Seeded row ordering exercises consumers that must not rely on table order.
    rng = np.random.default_rng(seed)
    produced: list[str] = []
    for name, frame in frames.items():
        filename = f"{name}_results.parquet"
        path = output_dir / filename
        if path.exists():
            raise FileExistsError(f"Refusing to overwrite existing fixture: {path}")
        order = rng.permutation(len(frame))
        frame.iloc[order].to_parquet(path, index=False, engine="pyarrow")
        produced.append(filename)
    return produced


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    names = build_fixtures(args.output_dir, args.seed)
    print(f"Built {len(names)} fixtures in {args.output_dir}")
    for name in names:
        print(name)


if __name__ == "__main__":
    main()
