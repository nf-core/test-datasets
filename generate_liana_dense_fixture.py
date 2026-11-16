"""Generate a dense synthetic fixture with complete LIANA consensus coverage."""

from __future__ import annotations

import argparse
import hashlib
import platform
from importlib.metadata import distribution, version
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from generate_de_raw_fixtures import require, validate_hdf5

SEED = 20261025
FILENAME = "liana_dense_consensus.h5ad"
RESOURCE_SHA256 = "9b4880f638198725177805d1cf5777b8d2f512444613d579158f1529c196ca2a"
PAIRS = (
    ("CXCL12", "CXCR4"),
    ("CCL5", "CCR5"),
    ("TGFB1", "TGFBR1_TGFBR2"),
    ("IL6", "IL6R_IL6ST"),
    ("MIF", "CD74_CXCR4"),
)


def consensus_resource() -> tuple[pd.DataFrame, list[str]]:
    require(version("liana") == "1.6.1", "The resource requires liana==1.6.1")
    path = distribution("liana").locate_file("liana/resource/omni_resource.csv")
    require(
        hashlib.sha256(path.read_bytes()).hexdigest() == RESOURCE_SHA256,
        "Bundled LIANA resource hash differs",
    )
    resource = pd.read_csv(path)
    resource = resource.loc[
        resource["resource"] == "consensus", ["source_genesymbol", "target_genesymbol"]
    ].rename(columns={"source_genesymbol": "ligand", "target_genesymbol": "receptor"})
    genes = sorted(
        {gene for entry in resource.to_numpy().ravel() for gene in entry.split("_")}
    )
    require(
        resource.shape == (4624, 2) and len(genes) == 1877,
        "Unexpected consensus resource",
    )
    pairs = set(resource.itertuples(index=False, name=None))
    require(all(pair in pairs for pair in PAIRS), "A documented interaction is absent")
    return resource, genes


def make_fixture(genes: list[str]) -> ad.AnnData:
    rng = np.random.default_rng(SEED)
    counts = rng.integers(1, 4, size=(40, len(genes)), dtype=np.int32)
    index = {gene: i for i, gene in enumerate(genes)}
    ligands = sorted({gene for ligand, _ in PAIRS for gene in ligand.split("_")})
    receptors = sorted({gene for _, receptor in PAIRS for gene in receptor.split("_")})
    for rows, expressed in ((slice(0, 20), ligands), (slice(20, 40), receptors)):
        columns = [index[gene] for gene in expressed]
        counts[rows, columns] = rng.integers(
            24, 28, size=(20, len(columns)), dtype=np.int32
        )
    obs = pd.DataFrame(
        {
            "group": pd.Categorical(
                ["sender"] * 20 + ["receiver"] * 20, categories=["sender", "receiver"]
            )
        },
        index=pd.Index([f"cell_{i:03d}" for i in range(40)], name="cell_id"),
    )
    var = pd.DataFrame(
        {"symbols": np.asarray(genes, dtype=object)},
        index=pd.Index(genes, name="gene_id"),
    )
    return ad.AnnData(X=counts, obs=obs, var=var)


def validate_fixture(data: ad.AnnData, genes: list[str]) -> None:
    require(data.shape == (40, 1877), "Unexpected fixture shape")
    require(
        isinstance(data.X, np.ndarray) and data.X.dtype == np.int32,
        "X must be dense int32",
    )
    require(
        np.isfinite(data.X).all() and (data.X > 0).all(),
        "All resource genes must be expressed",
    )
    require(data.obs_names.is_unique and data.var_names.is_unique, "Duplicate axes")
    require(
        data.var_names.tolist() == genes, "Resource component gene coverage differs"
    )
    require(data.var["symbols"].tolist() == genes, "Gene symbols differ")
    require(
        data.obs["group"].value_counts().to_dict() == {"sender": 20, "receiver": 20},
        "Groups must each contain 20 cells",
    )
    require(
        data.raw is None and not data.layers, "Unexpected alternative expression source"
    )
    sender = data.obs["group"] == "sender"
    for ligand, receptor in PAIRS:
        for component in ligand.split("_"):
            values = data.X[:, data.var_names.get_loc(component)]
            require(
                (values[sender] >= 24).all() and (values[~sender] <= 3).all(),
                f"Ligand {component} lacks sender enrichment",
            )
        for component in receptor.split("_"):
            values = data.X[:, data.var_names.get_loc(component)]
            require(
                (values[~sender] >= 24).all() and (values[sender] <= 3).all(),
                f"Receptor {component} lacks receiver enrichment",
            )


def validate_liana(data: ad.AnnData) -> None:
    import liana as li
    import scanpy as sc
    from threadpoolctl import threadpool_limits

    working = data.copy()
    with threadpool_limits(limits=1):
        sc.pp.log1p(working)
        li.mt.rank_aggregate(
            working,
            "group",
            resource_name="consensus",
            use_raw=False,
            n_jobs=1,
            n_perms=5,
            seed=SEED,
            verbose=True,
        )
    result = working.uns["liana_res"]
    require(not result.empty, "LIANA returned no interactions")
    require(
        len(result) == 4 * 4624,
        "Expected every consensus pair for all four group combinations",
    )
    print(f"LIANA result shape: {result.shape}; columns: {list(result.columns)}")
    for column in ("magnitude_rank", "specificity_rank"):
        require(
            np.isfinite(result[column]).all() and result[column].between(0, 1).all(),
            f"Invalid {column}",
        )
    directed = result.loc[
        (result["source"] == "sender") & (result["target"] == "receiver")
    ]
    for ligand, receptor in PAIRS:
        rows = directed.loc[
            (directed["ligand_complex"] == ligand)
            & (directed["receptor_complex"] == receptor)
        ]
        require(
            len(rows) == 1,
            f"Expected one sender-to-receiver {ligand}/{receptor} interaction",
        )
        require(
            (rows["lr_means"] > 3).all() and (rows["lr_logfc"] > 1).all(),
            f"Expected enriched expression and positive fold change for {ligand}/{receptor}",
        )
        print(
            rows[
                [
                    "source",
                    "target",
                    "ligand_complex",
                    "receptor_complex",
                    "lr_means",
                    "lr_logfc",
                    "magnitude_rank",
                    "specificity_rank",
                ]
            ].to_string(index=False)
        )
    print(
        f"Validated {len(PAIRS)} documented sender-to-receiver interactions; scanpy=={version('scanpy')}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--validate-liana",
        action="store_true",
        help="Run real LIANA 1.6.1 rank aggregation with five permutations",
    )
    args = parser.parse_args()
    path = args.output_dir / FILENAME
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite {path}")
    _, genes = consensus_resource()
    data = make_fixture(genes)
    validate_fixture(data, genes)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    data.write_h5ad(path, compression="gzip", compression_opts=4)
    validate_hdf5(path)
    restored = ad.read_h5ad(path)
    validate_fixture(restored, genes)
    np.testing.assert_array_equal(data.X, restored.X)
    pd.testing.assert_frame_equal(data.obs, restored.obs)
    pd.testing.assert_frame_equal(data.var, restored.var)
    print(f"Validated dense {data.n_obs}-cell/{data.n_vars}-gene fixture; seed={SEED}")
    print(
        f"Consensus component coverage: {len(genes)}/{len(genes)}; missing proportion=0"
    )
    print(f"Resource SHA-256: {RESOURCE_SHA256}")
    print(f"python=={platform.python_version()}")
    for package in ("liana", "anndata", "numpy", "pandas", "scipy", "h5py"):
        print(f"{package}=={version(package)}")
    print(f"{path}: {path.stat().st_size} bytes")
    if args.validate_liana:
        validate_liana(restored)


if __name__ == "__main__":
    main()
