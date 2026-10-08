"""Generate small DE tables for the plotting regression scenarios."""

from pathlib import Path

import pandas as pd

outdir = Path(__file__).parent


def table(groups):
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
            for gene, fold_change, pvalue in [("GENE1", 2.0, 0.001), ("GENE2", -1.5, 0.01), ("GENE3", 0.2, 0.5)]
        ]
    )


two_groups = table(["B", "A"])
two_groups_one_panel = two_groups.copy()
two_groups_one_panel.loc[two_groups_one_panel["group"] == "A", ["pvalue", "padj"]] = float("nan")
three_groups = table(["C", "B", "A"])
sparse_groups = three_groups.copy()
sparse_groups.loc[sparse_groups["group"] == "C", ["pvalue", "padj"]] = float("nan")
one_panel = sparse_groups.copy()
one_panel.loc[one_panel["group"] == "A", ["pvalue", "padj"]] = float("nan")
missing_pvalues = table(["A"])
missing_pvalues[["pvalue", "padj"]] = float("nan")
raw_pvalues = table(["A"])
raw_pvalues["padj"] = float("nan")

for name, frame in {
    "two_groups": two_groups,
    "two_groups_one_panel": two_groups_one_panel,
    "three_groups": three_groups,
    "sparse_groups": sparse_groups,
    "one_panel": one_panel,
    "missing_pvalues": missing_pvalues,
    "raw_pvalues": raw_pvalues,
    "condition:leiden:11:wilcoxon": two_groups,
}.items():
    frame.to_parquet(outdir / f"{name}_results.parquet", index=False, engine="pyarrow")
