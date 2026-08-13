#!/usr/bin/env python3
"""Validate the semantic contract of the compact nf-core/gwas fixtures."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import statistics
import struct
from collections import Counter
from pathlib import Path
from urllib.parse import urlparse

from generate_gwas_fixtures import ANALYSIS_HEADER, FIXTURE_BASE_URL

GENOTYPE_PATTERN = re.compile(r"^[01]\|[01]$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--pheno", type=Path, required=True)
    parser.add_argument("--qcovar", type=Path, required=True)
    parser.add_argument("--catcovar", type=Path, required=True)
    parser.add_argument("--relational-dir", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=200)
    parser.add_argument("--chromosomes", default="1,2")
    parser.add_argument("--variants-per-chromosome", type=int, default=1100)
    parser.add_argument("--cases", type=int, default=87)
    return parser.parse_args()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def assert_bgzf(path: Path) -> None:
    header = path.read_bytes()[:18]
    require(len(header) == 18, f"{path}: truncated gzip header")
    require(header[:3] == b"\x1f\x8b\x08", f"{path}: not gzip deflate data")
    require(header[3] & 4, f"{path}: gzip extra field is absent")
    extra_length = struct.unpack("<H", header[10:12])[0]
    require(extra_length >= 6, f"{path}: BGZF extra field is too short")
    require(header[12:16] == b"BC\x02\x00", f"{path}: BGZF BC field is absent")
    with gzip.open(path, "rb") as handle:
        while handle.read(1024 * 1024):
            pass


def read_tsv(path: Path, expected_header: list[str]) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames == expected_header, f"{path}: unexpected header")
        rows = list(reader)
    require(all(row["FID"] == row["IID"] for row in rows), f"{path}: FID/IID mismatch")
    return rows


def read_csv(path: Path, expected_header: list[str]) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == expected_header, f"{path}: unexpected header")
        return list(reader)


def fixture_url_to_path(
    url: str, fixtures_dir: Path, generated_files: dict[str, Path]
) -> Path:
    prefix = f"{FIXTURE_BASE_URL}/"
    require(url.startswith(prefix), f"URL is not on the stable gwas branch: {url}")
    parsed = urlparse(url)
    require(
        not parsed.query and not parsed.fragment, f"URL has query or fragment: {url}"
    )
    relative_path = url.removeprefix(prefix)
    require(".." not in Path(relative_path).parts, f"URL escapes fixture root: {url}")
    local_path = generated_files.get(relative_path, fixtures_dir / relative_path)
    require(local_path.is_file(), f"URL target is absent from fixture bundle: {url}")
    return local_path


def validate_relational_fixtures(
    relational_dir: Path,
    variant_ids: set[str],
    generated_files: dict[str, Path],
) -> dict[str, object]:
    fixtures_dir = relational_dir.parent
    expected_files = {
        "analysis_manifest_association_only.csv",
        "analysis_manifest_binary.csv",
        "analysis_manifest_heritability_only.csv",
        "analysis_manifest_heterogeneous.csv",
        "analysis_manifest_quantitative.csv",
        "cohort_manifest.csv",
        "method_options_heterogeneous.json",
        "resources/gcta_grm_extract.txt",
        "resources/ldak_predictor_extract.txt",
        "resources/ldak_weights.txt",
    }
    actual_files = {
        path.relative_to(relational_dir).as_posix()
        for path in relational_dir.rglob("*")
        if path.is_file()
    }
    require(actual_files == expected_files, "unexpected relational fixture inventory")

    cohort_header = [
        "cohort_id",
        "genome_build",
        "ancestry",
        "pgen",
        "psam",
        "pvar",
        "bed",
        "bim",
        "fam",
        "vcf",
    ]
    cohort_rows = read_csv(relational_dir / "cohort_manifest.csv", cohort_header)
    require(len(cohort_rows) == 1, "cohort manifest must contain exactly one row")
    cohort = cohort_rows[0]
    require(cohort["cohort_id"] == "example_cohort", "unexpected cohort ID")
    require(cohort["genome_build"] == "GRCh37", "unexpected genome build")
    require(cohort["ancestry"] == "EUR", "unexpected ancestry")
    require(
        all(
            not cohort[column]
            for column in ("pgen", "psam", "pvar", "bed", "bim", "fam")
        ),
        "cohort manifest must select only the canonical VCF representation",
    )
    fixture_url_to_path(cohort["vcf"], fixtures_dir, generated_files)

    expected_analyses = {
        "analysis_manifest_quantitative.csv": {
            "example_quantitative": ("QT", "quantitative", {"plink2"}, set()),
        },
        "analysis_manifest_binary.csv": {
            "example_binary": ("BT", "binary", {"plink2"}, set()),
        },
        "analysis_manifest_association_only.csv": {
            "example_association": (
                "QT",
                "quantitative",
                {"plink2", "regenie", "gcta_fastgwa", "ldak_kvik"},
                set(),
            ),
        },
        "analysis_manifest_heritability_only.csv": {
            "example_heritability": (
                "BT",
                "binary",
                set(),
                {"gcta_greml", "gcta_greml_ldms", "ldak_reml", "ldak_he", "ldak_pcgc"},
            ),
        },
        "analysis_manifest_heterogeneous.csv": {
            "heterogeneous_qt": (
                "QT",
                "quantitative",
                {"plink2", "regenie", "gcta_fastgwa", "ldak_kvik"},
                {"gcta_greml", "gcta_greml_ldms", "ldak_reml", "ldak_he"},
            ),
            "heterogeneous_bt": (
                "BT",
                "binary",
                {"plink2", "gcta_fastgwa", "ldak_kvik"},
                {"gcta_greml", "gcta_greml_ldms", "ldak_reml", "ldak_he", "ldak_pcgc"},
            ),
        },
    }
    all_analysis_ids: set[str] = set()
    for filename, expected_rows in expected_analyses.items():
        rows = read_csv(relational_dir / filename, ANALYSIS_HEADER)
        require(
            {row["analysis_id"] for row in rows} == set(expected_rows),
            f"{filename}: unexpected analysis IDs",
        )
        for row in rows:
            analysis_id = row["analysis_id"]
            require(
                analysis_id not in all_analysis_ids,
                f"duplicate analysis ID {analysis_id}",
            )
            all_analysis_ids.add(analysis_id)
            trait_id, trait_type, associations, heritability = expected_rows[
                analysis_id
            ]
            require(
                row["cohort_id"] == cohort["cohort_id"],
                f"{analysis_id}: unknown cohort",
            )
            require(
                (row["trait_id"], row["trait_type"], row["phenotype_column"])
                == (trait_id, trait_type, trait_id),
                f"{analysis_id}: inconsistent trait contract",
            )
            fixture_url_to_path(row["phenotype"], fixtures_dir, generated_files)
            for column in ("quant_covariates", "cat_covariates"):
                if row[column]:
                    fixture_url_to_path(row[column], fixtures_dir, generated_files)
            require(
                set(filter(None, row["association_methods"].split(",")))
                == associations,
                f"{analysis_id}: unexpected association methods",
            )
            require(
                set(filter(None, row["heritability_methods"].split(",")))
                == heritability,
                f"{analysis_id}: unexpected heritability methods",
            )
            if trait_type == "binary":
                require(
                    (row["control_value"], row["case_value"]) == ("1", "2"),
                    f"{analysis_id}: invalid binary coding",
                )
            else:
                require(
                    not row["control_value"] and not row["case_value"],
                    f"{analysis_id}: quantitative trait has binary coding",
                )
            if "ldak_pcgc" in heritability:
                require(
                    row["population_prevalence"] == "0.1",
                    f"{analysis_id}: invalid prevalence",
                )

    with (relational_dir / "method_options_heterogeneous.json").open(
        encoding="utf-8"
    ) as handle:
        method_options = json.load(handle)
    require(
        set(method_options) == {"heterogeneous_qt", "heterogeneous_bt"},
        "method options do not join to the heterogeneous analyses",
    )
    require(
        all(set(options) == {"gcta", "ldak"} for options in method_options.values()),
        "method options have unexpected method keys",
    )
    for option_name in ("grm_extract",):
        fixture_url_to_path(
            method_options["heterogeneous_qt"]["gcta"][option_name],
            fixtures_dir,
            generated_files,
        )
    for option_name in ("weights", "predictor_extract"):
        fixture_url_to_path(
            method_options["heterogeneous_qt"]["ldak"][option_name],
            fixtures_dir,
            generated_files,
        )

    selected_ids: set[str] = set()
    for filename in ("gcta_grm_extract.txt", "ldak_predictor_extract.txt"):
        ids = (
            (relational_dir / "resources" / filename)
            .read_text(encoding="utf-8")
            .splitlines()
        )
        require(
            ids and len(ids) == len(set(ids)), f"{filename}: empty or duplicate IDs"
        )
        require(set(ids) <= variant_ids, f"{filename}: IDs are absent from compact VCF")
        selected_ids.update(ids)
    weight_lines = (
        (relational_dir / "resources" / "ldak_weights.txt")
        .read_text(encoding="utf-8")
        .splitlines()
    )
    weights = [line.split() for line in weight_lines]
    require(
        weights
        and all(len(fields) == 2 and float(fields[1]) > 0 for fields in weights),
        "ldak_weights.txt: invalid predictor weights",
    )
    require(
        {fields[0] for fields in weights} <= variant_ids,
        "ldak_weights.txt: IDs are absent from compact VCF",
    )
    selected_ids.update(fields[0] for fields in weights)

    return {
        "relational_files": len(actual_files),
        "cohorts": len(cohort_rows),
        "analyses": len(all_analysis_ids),
        "resource_variant_ids": sorted(selected_ids),
    }


def squared_correlation(first: list[int], second: list[int]) -> float:
    first_mean = sum(first) / len(first)
    second_mean = sum(second) / len(second)
    covariance = sum(
        (left - first_mean) * (right - second_mean)
        for left, right in zip(first, second, strict=True)
    )
    first_ss = sum((value - first_mean) ** 2 for value in first)
    second_ss = sum((value - second_mean) ** 2 for value in second)
    require(first_ss > 0 and second_ss > 0, "constant genotype vector")
    return covariance**2 / (first_ss * second_ss)


def parse_vcf(
    path: Path,
    expected_chromosomes: set[str],
) -> tuple[list[str], Counter[str], list[str], list[list[int]], list[int]]:
    sample_ids: list[str] = []
    chromosome_counts: Counter[str] = Counter()
    variant_ids: list[str] = []
    dosages: list[list[int]] = []
    positions: list[int] = []
    previous_key: tuple[int, int] | None = None

    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "#CHROM":
                require(
                    fields[1:9]
                    == ["POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"],
                    f"{path}: malformed VCF header",
                )
                sample_ids = fields[9:]
                continue

            require(bool(sample_ids), f"{path}: data before #CHROM header")
            require(len(fields) == 9 + len(sample_ids), f"{path}: malformed record")
            chromosome, position_text, variant_id, ref, alt, qual, filt, info, fmt = (
                fields[:9]
            )
            require(
                chromosome in expected_chromosomes,
                f"{path}: unexpected chromosome {chromosome}",
            )
            position = int(position_text)
            key = (int(chromosome), position)
            require(
                previous_key is None or key > previous_key,
                f"{path}: records are not sorted",
            )
            previous_key = key
            require(
                ref in "ACGT" and alt in "ACGT" and ref != alt,
                f"{path}: invalid alleles",
            )
            require("," not in alt, f"{path}: multiallelic record")
            require(
                (qual, filt, info, fmt) == (".", "PASS", ".", "GT"),
                f"{path}: record is not GT-only",
            )
            genotypes = fields[9:]
            require(
                all(GENOTYPE_PATTERN.fullmatch(value) for value in genotypes),
                f"{path}: invalid or missing genotype",
            )
            dosage = [int(value[0]) + int(value[2]) for value in genotypes]
            allele_count = sum(dosage)
            require(
                20 <= allele_count <= 2 * len(sample_ids) - 20,
                f"{path}: unusable variant {variant_id}",
            )

            chromosome_counts[chromosome] += 1
            variant_ids.append(variant_id)
            dosages.append(dosage)
            positions.append(position)

    require(len(variant_ids) == len(set(variant_ids)), f"{path}: duplicate variant IDs")
    require(
        all(len(variant_id) <= 8 for variant_id in variant_ids),
        f"{path}: variant ID is not short",
    )
    return sample_ids, chromosome_counts, variant_ids, dosages, positions


def matrix_rank(matrix: list[list[float]], tolerance: float = 1e-10) -> int:
    work = [row[:] for row in matrix]
    rows = len(work)
    columns = len(work[0])
    pivot_row = 0
    for column in range(columns):
        pivot = max(range(pivot_row, rows), key=lambda row: abs(work[row][column]))
        if abs(work[pivot][column]) <= tolerance:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        pivot_value = work[pivot_row][column]
        work[pivot_row] = [value / pivot_value for value in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row:
                continue
            factor = work[row][column]
            if abs(factor) > tolerance:
                work[row] = [
                    value - factor * pivot_value
                    for value, pivot_value in zip(
                        work[row], work[pivot_row], strict=True
                    )
                ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def validate() -> dict[str, object]:
    args = parse_args()
    chromosomes = args.chromosomes.split(",")
    assert_bgzf(args.vcf)
    samples, chromosome_counts, variant_ids, dosages, _ = parse_vcf(
        args.vcf, set(chromosomes)
    )

    expected_sample_ids = [f"S{index:03d}" for index in range(1, args.samples + 1)]
    require(samples == expected_sample_ids, "VCF sample IDs or order are unexpected")
    require(
        chromosome_counts
        == Counter(
            {chromosome: args.variants_per_chromosome for chromosome in chromosomes}
        ),
        "chromosome variant counts are unexpected",
    )
    require(
        len(variant_ids) == len(chromosomes) * args.variants_per_chromosome,
        "total variant count is unexpected",
    )

    pheno = read_tsv(args.pheno, ["FID", "IID", "QT", "BT"])
    qcovar = read_tsv(args.qcovar, ["FID", "IID", "AGE", "Q1", "Q2", "Q3"])
    catcovar = read_tsv(args.catcovar, ["FID", "IID", "SEX"])
    for path, rows in (
        (args.pheno, pheno),
        (args.qcovar, qcovar),
        (args.catcovar, catcovar),
    ):
        require(
            [row["IID"] for row in rows] == samples,
            f"{path}: sample IDs or order differ from VCF",
        )

    qt = [float(row["QT"]) for row in pheno]
    require(statistics.variance(qt) > 0.1, "quantitative phenotype is degenerate")
    binary_counts = Counter(int(row["BT"]) for row in pheno)
    require(
        binary_counts == Counter({1: args.samples - args.cases, 2: args.cases}),
        "binary phenotype split is unexpected",
    )

    quantitative_matrix = [
        [1.0, *[float(row[column]) for column in ("AGE", "Q1", "Q2", "Q3")]]
        for row in qcovar
    ]
    rank = matrix_rank(quantitative_matrix)
    require(rank == 5, f"quantitative covariate design rank is {rank}, expected 5")
    categories = Counter(row["SEX"] for row in catcovar)
    require(
        categories == Counter({"1": args.samples // 2, "2": args.samples // 2}),
        "categorical covariate is unbalanced or degenerate",
    )
    combined_matrix = [
        [*quantitative_row, 1.0 if catcovar[index]["SEX"] == "2" else 0.0]
        for index, quantitative_row in enumerate(quantitative_matrix)
    ]
    combined_rank = matrix_rank(combined_matrix)
    require(
        combined_rank == 6,
        f"combined covariate design rank is {combined_rank}, expected 6",
    )

    within_block_r2: list[float] = []
    between_block_r2: list[float] = []
    block_size = 11
    for chromosome_index in range(len(chromosomes)):
        chromosome_start = chromosome_index * args.variants_per_chromosome
        for block_start in range(
            chromosome_start,
            chromosome_start + args.variants_per_chromosome,
            block_size,
        ):
            if block_start + 1 < chromosome_start + args.variants_per_chromosome:
                within_block_r2.append(
                    squared_correlation(dosages[block_start], dosages[block_start + 1])
                )
            next_block = block_start + block_size
            if next_block < chromosome_start + args.variants_per_chromosome:
                between_block_r2.append(
                    squared_correlation(dosages[block_start], dosages[next_block])
                )

    median_within = statistics.median(within_block_r2)
    median_between = statistics.median(between_block_r2)
    require(
        0.05 < median_within < 0.8,
        f"within-block LD is degenerate ({median_within:.4f})",
    )
    require(
        median_within > 3 * median_between,
        "within-block LD is not distinguishable from background",
    )

    summary = {
        "samples": len(samples),
        "chromosome_variants": dict(chromosome_counts),
        "variants": len(variant_ids),
        "controls": binary_counts[1],
        "cases": binary_counts[2],
        "qt_variance": round(statistics.variance(qt), 6),
        "quantitative_covariate_rank": rank,
        "combined_covariate_rank": combined_rank,
        "categorical_counts": dict(categories),
        "median_within_block_r2": round(median_within, 6),
        "median_between_block_r2": round(median_between, 6),
    }
    generated_files = {
        f"genotypes/{args.vcf.name}": args.vcf,
        f"pheno_cov/{args.pheno.name}": args.pheno,
        f"pheno_cov/{args.qcovar.name}": args.qcovar,
        f"pheno_cov/{args.catcovar.name}": args.catcovar,
    }
    summary.update(
        validate_relational_fixtures(
            args.relational_dir, set(variant_ids), generated_files
        )
    )
    return summary


def main() -> None:
    summary = validate()
    for key, value in summary.items():
        print(f"{key}\t{value}")


if __name__ == "__main__":
    main()
