# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = ["pandas==3.0.0", "pyarrow==23.0.1"]
# ///
"""Materialise the doublet regression tables introduced by scdownstream PRs 393 and 395."""

import hashlib
import json
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq


def generate(root):
    destination = root / "intermediates/doublet_detection/regressions"
    destination.mkdir(parents=True, exist_ok=True)
    cells = [f"cell_{i:03d}" for i in range(128)]
    tables = {}

    def calls(name, positives, order=None, ids=None):
        order = list(range(127, -1, -1)) if order is None else order
        ids = cells if ids is None else ids
        tables[name] = pd.DataFrame(
            {"arbitrary_prefix": [i in positives for i in order]},
            index=[ids[i] for i in order],
        )

    calls("doubletdetection", {0, 1, 2})
    calls("scdblfinder", {1, 3})
    calls("scrublet", {0, 3})
    calls("solo", {4})
    calls("solo_single", {1})
    calls("all_true", set(range(128)))
    calls("zero", set())
    calls("zero_solo", set())
    calls("small", {0})
    calls("small_solo", {0})
    tables["existing"] = pd.DataFrame({"solo_class": ["singlet"]}, index=cells[:1])
    tables["duplicate"] = pd.DataFrame(
        {"call": [True, False]}, index=[cells[0], cells[0]]
    )
    tables["missing"] = pd.DataFrame({"call": False}, index=cells[1:])
    tables["extra"] = pd.DataFrame({"call": False}, index=cells + ["cell_128"])
    tables["columns"] = pd.DataFrame({"call": False, "score": 9999}, index=cells)
    tables["null"] = pd.DataFrame(
        {"call": pd.array([False] * 128, dtype="boolean")}, index=cells
    )
    tables["null"].loc[cells[5], "call"] = pd.NA
    tables["invalid"] = pd.DataFrame(
        {"call": ["singlet" if i == 5 else "false" for i in range(128)]}, index=cells
    )
    tables["numeric"] = pd.DataFrame({"call": [i % 2 for i in range(128)]}, index=cells)

    order = list(range(64, 128)) + list(range(64))
    tables["doubletdetection_scores"] = pd.DataFrame(
        {"doubletdetection_score": [1000 + i for i in order]},
        index=[cells[i] for i in order],
    )
    tables["scdblfinder_scores"] = pd.DataFrame(
        {
            "scDblFinder.score": [2000 + i for i in order],
            "scDblFinder.class": [
                "doublet" if i in {1, 3} else "singlet" for i in order
            ],
            "scDblFinder.cxds_score": [i / 2 for i in order],
            "scDblFinder.mostLikelyOrigin": "A+B",
        },
        index=[cells[i] for i in order],
    )
    tables["scrublet_scores"] = pd.DataFrame(
        {"scrublet_score": [None if i == 2 else 3000 + i for i in order]},
        index=[cells[i] for i in order],
    )
    tables["scrublet_scores_full"] = pd.DataFrame(
        {"scrublet_score": [3000 + i for i in range(127, -1, -1)]},
        index=list(reversed(cells)),
    )

    for sample, offset, sc_calls, scr_calls in (
        ("sample_a", 1000, {1, 2, 7}, {0, 2, 7}),
        ("sample_b", 2000, {4, 5, 10}, {3, 5, 10}),
    ):
        ids = [f"{sample}_{cell}" for cell in cells]
        calls(f"{sample}.scdblfinder", sc_calls, ids=ids)
        calls(
            f"{sample}.scrublet",
            scr_calls,
            list(range(1, 128, 2)) + list(range(0, 128, 2)),
            ids,
        )
        order = list(range(127, -1, -1))
        tables[f"{sample}.scores"] = pd.DataFrame(
            {
                "scDblFinder.score": [offset + i + 0.5 for i in order],
                "scDblFinder.class": [
                    "doublet" if i in sc_calls else "singlet" for i in order
                ],
                "scrublet_score": [offset + i + 0.25 for i in order],
            },
            index=list(reversed(ids)),
        )

    manifest = {
        "generator": Path(__file__).name,
        "source": "scdownstream regression test literals from PRs 393 and 395",
        "tables": {},
    }
    for name, frame in sorted(tables.items()):
        table = frame.copy()
        table.insert(0, "__index__", frame.index)
        path = destination / f"{name}.parquet"
        table.to_parquet(path, engine="pyarrow", index=False, compression="snappy")
        pd.testing.assert_frame_equal(
            pd.read_parquet(path).set_index("__index__").rename_axis(None), frame
        )
        metadata = pq.read_metadata(path)
        assert all(
            metadata.row_group(0).column(i).compression == "SNAPPY"
            for i in range(metadata.num_columns)
        )
        manifest["tables"][path.relative_to(root).as_posix()] = {
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "rows": len(frame),
            "columns": list(table.columns),
        }
    (destination / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Verified {len(tables)} additive Parquet fixtures in {destination}")


if __name__ == "__main__":
    generate(Path(__file__).resolve().parent)
