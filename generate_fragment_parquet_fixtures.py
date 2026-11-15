# /// script
# requires-python = ">=3.12,<3.13"
# dependencies = ["numpy==2.3.5", "pandas==3.0.0", "pyarrow==23.0.1"]
# ///
"""Generate and verify supplemental scdownstream fragment fixtures."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

FRAGMENTS = Path("intermediates/fragments")
REGRESSIONS = FRAGMENTS / "regressions"
MANIFEST = FRAGMENTS / "parquet_manifest.json"
CELL = "SRR28679756::CCTTTCTTCTGCTGCT"
SECOND_CELL = "SRR28679756::ATCATCTTCTGCTGTC"
TESTS = {
    "extend": "modules/local/adata/extend/tests/main.nf.test",
    "entropy": "modules/local/adata/entropy/tests/main.nf.test",
    "splitcol": "modules/local/adata/splitcol/tests/labels.nf.test",
    "annotations": "modules/local/custom/aggregate_per_cell_annotations/tests/main.nf.test",
    "doublets": "modules/local/custom/doubletremoval/tests/main.nf.test",
    "doublet_plots": "modules/local/custom/doubletremoval/tests/plot_edge_cases.nf.test",
    "edgepython": "modules/local/edgepython/differential/tests/design.nf.test",
    "liana": "modules/local/liana/rankaggregate/tests/dense.nf.test",
    "scanorama": "modules/local/scanorama/integrate/tests/main.nf.test",
    "batch_validation": "modules/local/scanorama/integrate/tests/batch_validation.nf.test",
    "panels": "modules/local/scanpy/filter/tests/small_panels.nf.test",
    "panel_plots": "modules/local/scanpy/filter/tests/small_panels_plot.nf.test",
    "finalize": "subworkflows/local/finalize/tests/main.nf.test",
}


def csv_frame(payload: str) -> pd.DataFrame:
    return pd.read_csv(io.StringIO(payload), index_col=0)


def physical_frame(frame: pd.DataFrame, slot: str) -> pd.DataFrame:
    result = frame.copy()
    if slot != "uns":
        if not result.index.is_unique or result.index.hasnans:
            raise ValueError("Fragment identifiers must be unique and non-null")
        if not all(isinstance(value, str) for value in result.index):
            raise ValueError("Fragment identifiers must be strings")
    result.columns = result.columns.map(str)
    if "__index__" in result.columns:
        raise ValueError("Reserved __index__ column already exists")
    result.index.name = "__index__"
    return result.reset_index()


def parquet_bytes(frame: pd.DataFrame) -> bytes:
    buffer = io.BytesIO()
    frame.to_parquet(buffer, index=False, engine="pyarrow", compression="snappy")
    return buffer.getvalue()


def fixtures(root: Path) -> tuple[dict, dict, list[str]]:
    tables = {}
    text_files = {}
    retained = []
    required = [
        FRAGMENTS / f"{stem}.pkl"
        for stem in (
            "X_scvi",
            "X_scvi-global_umap",
            "scvi-global-1.0_entropy",
            "scvi-global-1.0_leiden",
            "scvi-global-1.0_paga",
        )
    ]
    required.extend(
        Path("intermediates/doublet_detection") / f"test_{name}.pkl"
        for name in ("doubletdetection", "scrublet", "solo")
    )
    missing = [str(path) for path in required if not (root / path).is_file()]
    if missing:
        raise ValueError("Missing legacy fixtures: " + ", ".join(missing))

    def add(name, frame, slot, consumers, source, malformed=None):
        relative = REGRESSIONS / f"{name}.parquet"
        tables[relative] = (frame, slot, consumers, source, malformed)

    for directory in (FRAGMENTS, Path("intermediates/doublet_detection")):
        for source in sorted((root / directory).glob("*.pkl")):
            value = pd.read_pickle(source)
            relative = source.relative_to(root)
            if not isinstance(value, pd.DataFrame):
                retained.append(relative.as_posix())
                continue
            slot = (
                "pred"
                if directory.name == "doublet_detection"
                else ("obsm" if source.stem.startswith("X_") else "obs")
            )
            consumers = ["doublets"] if slot == "pred" else ["extend", "finalize"]
            tables[relative.with_suffix(".parquet")] = (
                physical_frame(value, slot),
                slot,
                consumers,
                relative.as_posix(),
                None,
            )

    literals = {
        "a_obs": (
            f"cell,label,n_counts,score,flag,version.1\n{CELL},T_cell,100,1.0,True,first\n{SECOND_CELL},,,2.0,,\nnot_a_cell,ignored,999,9,True,ignored\n",
            "obs",
            "extend",
        ),
        "b_obs": (
            f"cell,label,n_counts,score,flag,version.1\n{CELL},T_cell,100.0,1,True,first\n{SECOND_CELL},B_cell,,2,False,second\n",
            "obs",
            "extend",
        ),
        "a_var": ("gene,selected,score\nFAM177A1P1,True,1.0\n", "var", "extend"),
        "b_var": (
            "gene,selected,score\nFAM177A1P1,True,1\nPOMP,False,2.5\n",
            "var",
            "extend",
        ),
        "a_conflict": (f"cell,label\n{CELL},T_cell\n", "obs", "extend"),
        "b_conflict": (f"cell,label\n{CELL},B_cell\n", "obs", "extend"),
        "a_score": ("gene,score\nFAM177A1P1,1\n", "var", "extend"),
        "b_score": ("gene,score\nFAM177A1P1,2\n", "var", "extend"),
        "annotations": (
            f"cell,annotation:test:per_cell\n{CELL},T_cell\n{SECOND_CELL},T_cell\n",
            "obs",
            "annotations",
        ),
        "missing_batch": (
            f"cell,regression_batch\n{CELL},batch_a\n",
            "obs",
            "batch_validation",
        ),
        "empty_batch": ("cell,empty_batch\nnot_a_cell,\n", "obs", "entropy"),
    }
    for name, (payload, slot, consumer) in literals.items():
        add(
            name,
            physical_frame(csv_frame(payload), slot),
            slot,
            [consumer],
            TESTS[consumer],
        )

    add(
        "duplicate_index",
        pd.DataFrame({"__index__": [CELL, CELL], "label": ["T_cell", "T_cell"]}),
        "obs",
        ["extend"],
        TESTS["extend"],
        "duplicate identifiers",
    )
    add(
        "missing_index",
        pd.DataFrame({"label": ["T_cell"]}),
        "obs",
        ["extend"],
        "Parquet reader contract",
        "missing __index__ column",
    )
    add(
        "null_index",
        pd.DataFrame(
            {
                "__index__": pd.Series([CELL, None], dtype="str"),
                "label": ["T_cell", "B_cell"],
            }
        ),
        "obs",
        ["extend"],
        "Parquet reader contract",
        "null identifier",
    )
    text_files[REGRESSIONS / "duplicate_columns.csv"] = (
        f"cell,label,label\n{CELL},T_cell,B_cell\n".encode(),
        ["extend"],
        TESTS["extend"],
    )

    cells = [f"cell_{i:03d}" for i in range(128)]
    labels = pd.DataFrame(
        {
            "safe": [
                ["A B", "a_b", "T:CD4+", "B-cell.(naive)"][i % 4] for i in range(128)
            ],
            "collision": [
                "safe" if i == 0 else "A B" if i % 2 else "A_B" for i in range(128)
            ],
            "path": ["safe"] + ["bad/name"] * 127,
            "traversal": ["safe"] + [".."] * 127,
            "staged_input": ["labels" if i % 2 else "safe" for i in range(128)],
        },
        index=cells,
    )
    add("labels", physical_frame(labels, "obs"), "obs", ["splitcol"], TESTS["splitcol"])

    for n_rows, name, consumer in (
        (128, "panels", "panels"),
        (32, "panels_plot", "panel_plots"),
    ):
        frame = pd.DataFrame(
            {f"panel_{n}": [i < n for i in range(n_rows)] for n in (1, 10, 19, 20, 32)},
            index=[f"GENE{i:03d}" for i in range(n_rows)],
        )
        add(name, physical_frame(frame, "var"), "var", [consumer], TESTS[consumer])
    frame = pd.DataFrame({"pct_counts_in_top_20_genes": [17.0] * 128}, index=cells)
    add("previous_qc", physical_frame(frame, "obs"), "obs", ["panels"], TESTS["panels"])

    matrix = pd.DataFrame(
        [[1 + (i + j) % 7 for j in range(32)] for i in range(8)],
        index=[f"cell_{i}" for i in range(8)],
        columns=[f"GENE{j:03d}" for j in range(32)],
    )
    matrix.index.name = "cell"
    add(
        "positive_counts",
        physical_frame(matrix, "counts"),
        "counts",
        ["panel_plots"],
        TESTS["panel_plots"],
    )
    text_files[REGRESSIONS / "positive_counts.csv"] = (
        matrix.to_csv(lineterminator="\n").encode(),
        ["panel_plots"],
        TESTS["panel_plots"],
    )

    for name, small in (("zero_predictions", False), ("small_predictions", True)):
        frame = pd.DataFrame(
            {
                column: [small and i == 0 for i in range(128)]
                for column in ("first", "second")
            },
            index=cells,
        )
        add(
            name,
            physical_frame(frame, "pred"),
            "pred",
            ["doublet_plots"],
            TESTS["doublet_plots"],
        )

    donor_rows = [
        (f"{donor}__{label}__{condition}", f"{donor}_{condition}")
        for donor in ("D01", "D02", "D03")
        for label in ("TypeA", "TypeB")
        for condition in ("Healthy", "Disease")
    ]
    frame = pd.DataFrame(
        {"unpaired_donor": [value for _, value in donor_rows]},
        index=[key for key, _ in donor_rows],
    )
    add(
        "unpaired_obs",
        physical_frame(frame, "obs"),
        "obs",
        ["edgepython"],
        TESTS["edgepython"],
    )

    csv_source = "pbmc/SRR28679756_filtered_matrix.csv"
    csv_cells = pd.read_csv(root / csv_source, usecols=[0]).iloc[:, 0].tolist()
    if len(csv_cells) != 64:
        raise ValueError(f"Expected 64 dense CSV cells, found {len(csv_cells)}")
    for mode in ("interleaved", "contiguous"):
        batches = [
            "z_batch"
            if (i % 2 == 0 if mode == "interleaved" else i < 32)
            else "a_batch"
            for i in range(len(csv_cells))
        ]
        frame = pd.DataFrame({"regression_batch": batches}, index=csv_cells)
        add(
            mode,
            physical_frame(frame, "obs"),
            "obs",
            ["scanorama"],
            [TESTS["scanorama"], csv_source],
        )
    frame = pd.DataFrame(
        {
            "interaction_group": [
                "sender" if i % 2 == 0 else "receiver" for i in range(len(csv_cells))
            ]
        },
        index=csv_cells,
    )
    add(
        "groups",
        physical_frame(frame, "obs"),
        "obs",
        ["liana"],
        [TESTS["liana"], csv_source],
    )

    large_integer = pd.DataFrame(
        {"large_integer": pd.array([9007199254740993, None], dtype="Int64")},
        index=[CELL, SECOND_CELL],
    )
    add(
        "large_integer",
        physical_frame(large_integer, "obs"),
        "obs",
        ["extend"],
        "Synthetic nullable integer precision contract",
    )

    uns = pd.DataFrame(
        {
            "source": ["sender", "receiver"],
            "target": ["receiver", "sender"],
            "ligand_complex": ["CXCL12", "CCL5"],
            "receptor_complex": ["CXCR4", "CCR5"],
            "magnitude_rank": [0.125, 0.5],
            "specificity_rank": [0.25, np.nan],
        },
        index=pd.Index([7, 3], dtype="int64"),
    )
    tables[REGRESSIONS / "uns_table.parquet"] = (
        physical_frame(uns, "uns"),
        "uns",
        ["extend", "finalize"],
        "Synthetic LIANA-like indexed table",
        None,
    )
    for slot, ids in (
        ("obs", [CELL, SECOND_CELL, "not_a_cell"]),
        ("var", ["FAM177A1P1", "POMP", "not_a_gene"]),
    ):
        frame = pd.DataFrame(
            {
                "ordered_label": pd.Categorical(
                    ["high", None, "low"],
                    categories=pd.Index(
                        ["low", "medium", "high", "unused"], dtype=object
                    ),
                    ordered=True,
                ),
                "unordered_label": pd.Categorical(
                    ["B", "A", None],
                    categories=pd.Index(["B", "A", "unused"], dtype=object),
                ),
                "nullable_integer": pd.array([1, None, 3], dtype="Int64"),
                "nullable_boolean": pd.array([True, None, False], dtype="boolean"),
                "nullable_string": pd.array(["001", None, "NA"], dtype="string"),
                "nullable_float": [1.5, np.nan, 0.0],
            },
            index=ids,
        )
        add(
            f"categorical_{slot}",
            physical_frame(frame, slot),
            slot,
            ["extend"],
            "Synthetic reader round-trip contract",
        )
    return tables, text_files, retained


def validate(
    frame: pd.DataFrame, payload: bytes, slot: str, malformed: str | None
) -> dict:
    restored = pd.read_parquet(io.BytesIO(payload), engine="pyarrow")
    # pandas 3 infers string category indexes when restoring object-backed strings.
    pd.testing.assert_frame_equal(
        restored, frame, check_exact=True, check_dtype=False, check_categorical=False
    )
    parquet = pq.ParquetFile(io.BytesIO(payload))
    schema = parquet.schema_arrow
    if not schema.equals(
        pa.Table.from_pandas(frame, preserve_index=False).schema, check_metadata=True
    ):
        raise ValueError("Arrow schema or pandas metadata differs from the source")
    for column in frame:
        if isinstance(frame[column].dtype, pd.CategoricalDtype):
            pd.testing.assert_index_equal(
                restored[column].cat.categories,
                frame[column].cat.categories,
                exact=False,
            )
            if restored[column].cat.ordered != frame[column].cat.ordered:
                raise ValueError("Categorical ordering differs")
            np.testing.assert_array_equal(
                restored[column].cat.codes, frame[column].cat.codes
            )
        elif (
            frame[column].dtype != object
            and restored[column].dtype != frame[column].dtype
        ):
            raise ValueError(f"Pandas dtype differs: {column}")
    metadata = json.loads(schema.metadata[b"pandas"])
    if metadata["index_columns"]:
        raise ValueError("Parquet must not contain a hidden pandas index")
    if list(schema.names) != list(frame.columns) or not schema.names:
        raise ValueError("Physical column order differs")
    if malformed is None:
        ids = restored["__index__"]
        if ids.isna().any() or not ids.is_unique:
            raise ValueError("Identifiers must be unique and non-null")
        if slot != "uns" and not all(isinstance(value, str) for value in ids):
            raise ValueError("Identifiers must be strings")
        if slot == "obsm" and list(frame.columns[1:]) != [
            str(i) for i in range(len(frame.columns) - 1)
        ]:
            raise ValueError("Embedding columns must be consecutive string integers")
    if any(
        parquet.metadata.row_group(i).column(j).compression != "SNAPPY"
        for i in range(parquet.metadata.num_row_groups)
        for j in range(len(schema))
    ):
        raise ValueError("Every column must use Snappy compression")
    fields = []
    for field in schema:
        entry = {
            "name": field.name,
            "arrow_type": str(field.type),
            "pandas_dtype": str(frame[field.name].dtype),
            "null_count": int(frame[field.name].isna().sum()),
        }
        if isinstance(frame[field.name].dtype, pd.CategoricalDtype):
            entry.update(
                categories=frame[field.name].cat.categories.tolist(),
                ordered=frame[field.name].cat.ordered,
            )
        fields.append(entry)
    return {"rows": len(frame), "fields": fields}


def build(root: Path, check: bool) -> None:
    tables, text_files, retained = fixtures(root)
    generated = {}
    manifest = {
        "schema_version": 1,
        "generator": "generate_fragment_parquet_fixtures.py",
        "source_pipeline_revision": "a219db8f144eb1f8be17104b8b4b733539396298",
        "runtime": {
            "python": "3.12",
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "pyarrow": pa.__version__,
        },
        "encoding": {
            "index_column": "__index__",
            "pandas_index": False,
            "engine": "pyarrow",
            "compression": "snappy",
        },
        "consumer_tests": TESTS,
        "retained_nested_pickles": retained,
        "source_sha256": {
            path.relative_to(root).as_posix(): hashlib.sha256(
                path.read_bytes()
            ).hexdigest()
            for path in sorted(
                [
                    *(root / FRAGMENTS).glob("*.pkl"),
                    *(root / "intermediates/doublet_detection").glob("*.pkl"),
                    root / "pbmc/SRR28679756_filtered_matrix.csv",
                ]
            )
        },
        "files": {},
    }
    for relative, (frame, slot, consumers, source, malformed) in sorted(tables.items()):
        payload = parquet_bytes(frame)
        if payload != parquet_bytes(frame):
            raise ValueError(f"Non-deterministic Parquet serialisation: {relative}")
        entry = validate(frame, payload, slot, malformed)
        entry.update(slot=slot, consumers=consumers, source=source)
        if malformed:
            entry["intentional_error"] = malformed
        generated[relative] = payload
        manifest["files"][relative.as_posix()] = entry
    for relative, (payload, consumers, source) in sorted(text_files.items()):
        generated[relative] = payload
        manifest["files"][relative.as_posix()] = {
            "format": "csv",
            "consumers": consumers,
            "source": source,
        }
    for relative, payload in generated.items():
        if len(payload) >= 1_000_000:
            raise ValueError(f"Fixture exceeds size limit: {relative}")
        manifest["files"][relative.as_posix()].update(
            bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest()
        )
    generated[MANIFEST] = (
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    ).encode()
    failures = []
    for relative, payload in generated.items():
        path = root / relative
        if check:
            if not path.is_file() or path.read_bytes() != payload:
                failures.append(relative.as_posix())
            elif relative.suffix == ".parquet":
                frame, slot, _, _, malformed = tables[relative]
                validate(frame, path.read_bytes(), slot, malformed)
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(payload)
    if failures:
        raise ValueError(
            "Missing or changed generated fixtures: " + ", ".join(failures)
        )
    print(
        f"{'Checked' if check else 'Generated'} {len(tables)} Parquet fixtures, {len(text_files)} CSV fixtures and {MANIFEST}"
    )
    print(
        f"Largest fixture: {max(len(payload) for path, payload in generated.items() if path != MANIFEST):,} bytes; byte determinism and value/schema round trips passed"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset-root",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Existing scdownstream dataset checkout containing the legacy fixtures",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Validate committed bytes and schemas without writing files",
    )
    args = parser.parse_args()
    build(args.dataset_root.resolve(), args.check)


if __name__ == "__main__":
    main()
