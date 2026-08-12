#!/usr/bin/env python3
"""Generate the compact, deterministic nf-core/gwas test fixtures."""

from __future__ import annotations

import argparse
import csv
import math
import random
from pathlib import Path

ALLELE_PAIRS = (("A", "C"), ("G", "T"), ("C", "G"), ("T", "A"))
BLOCK_SIZE = 11


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", default="example")
    parser.add_argument("--chromosomes", default="1,2")
    parser.add_argument("--samples", type=int, default=200)
    parser.add_argument("--variants-per-chromosome", type=int, default=1100)
    parser.add_argument("--cases", type=int, default=87)
    parser.add_argument("--seed", type=int, default=20260812)
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    return parser.parse_args()


def standardise(values: list[float]) -> list[float]:
    mean = sum(values) / len(values)
    variance = sum((value - mean) ** 2 for value in values) / (len(values) - 1)
    if variance == 0:
        raise ValueError("cannot standardise a constant vector")
    scale = math.sqrt(variance)
    return [(value - mean) / scale for value in values]


def normal_draws(rng: random.Random, count: int) -> list[float]:
    """Return fixed Box-Muller draws using only Random.random()."""
    draws: list[float] = []
    while len(draws) < count:
        u1 = max(rng.random(), 1e-15)
        u2 = rng.random()
        radius = math.sqrt(-2.0 * math.log(u1))
        angle = 2.0 * math.pi * u2
        draws.extend((radius * math.cos(angle), radius * math.sin(angle)))
    return draws[:count]


def write_tsv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def build_genotypes(
    rng: random.Random,
    chromosomes: list[int],
    samples: list[str],
    variants_per_chromosome: int,
) -> tuple[list[tuple[int, int, str, str, str, list[str]]], list[list[int]]]:
    records: list[tuple[int, int, str, str, str, list[str]]] = []
    dosages: list[list[int]] = []
    haplotype_count = 2 * len(samples)

    for chromosome in chromosomes:
        block_count = math.ceil(variants_per_chromosome / BLOCK_SIZE)
        for block_index in range(block_count):
            frequency_step = (block_index * 17 + chromosome * 13) % 100
            frequency = 0.10 + 0.35 * frequency_step / 99.0
            founder = [
                1 if rng.random() < frequency else 0 for _ in range(haplotype_count)
            ]

            for offset in range(BLOCK_SIZE):
                variant_number = block_index * BLOCK_SIZE + offset + 1
                if variant_number > variants_per_chromosome:
                    break
                if offset == 0:
                    alleles = founder
                else:
                    alleles = [
                        (
                            founder_allele
                            if rng.random() < 0.72
                            else (1 if rng.random() < frequency else 0)
                        )
                        for founder_allele in founder
                    ]

                ref, alt = ALLELE_PAIRS[
                    (variant_number + chromosome) % len(ALLELE_PAIRS)
                ]
                variant_id = f"v{chromosome}_{variant_number:04d}"
                position = 1_000 + variant_number * 100
                genotypes = [
                    f"{alleles[2 * index]}|{alleles[2 * index + 1]}"
                    for index in range(len(samples))
                ]
                variant_dosages = [
                    alleles[2 * index] + alleles[2 * index + 1]
                    for index in range(len(samples))
                ]
                records.append((chromosome, position, variant_id, ref, alt, genotypes))
                dosages.append(variant_dosages)

    return records, dosages


def write_vcf(
    path: Path,
    samples: list[str],
    records: list[tuple[int, int, str, str, str, list[str]]],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("##fileformat=VCFv4.3\n")
        handle.write("##source=nf-core-test-datasets-compact-gwas-fixture\n")
        handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        for chromosome in sorted({record[0] for record in records}):
            handle.write(f"##contig=<ID={chromosome}>\n")
        header = [
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
            *samples,
        ]
        handle.write("\t".join(header) + "\n")
        for chromosome, position, variant_id, ref, alt, genotypes in records:
            fields = [
                str(chromosome),
                str(position),
                variant_id,
                ref,
                alt,
                ".",
                "PASS",
                ".",
                "GT",
                *genotypes,
            ]
            handle.write("\t".join(fields) + "\n")


def polygenic_score(
    dosages: list[list[int]], start: int, step: int, count: int
) -> list[float]:
    sample_count = len(dosages[0])
    score = [0.0] * sample_count
    variant_indices = [(start + step * index) % len(dosages) for index in range(count)]
    for effect_index, variant_index in enumerate(variant_indices):
        magnitude = 0.15 + 0.05 * (effect_index % 7)
        effect = magnitude if effect_index % 2 == 0 else -magnitude
        for sample_index, dosage in enumerate(dosages[variant_index]):
            score[sample_index] += effect * dosage
    return standardise(score)


def write_sidecars(
    output_dir: Path,
    prefix: str,
    samples: list[str],
    dosages: list[list[int]],
    cases: int,
    rng: random.Random,
) -> None:
    sample_count = len(samples)
    if not 0 < cases < sample_count:
        raise ValueError("case count must be between zero and the sample count")

    age = [30 + (index * 7) % 40 for index in range(sample_count)]
    q1 = standardise([float(index) for index in range(sample_count)])
    q2 = [
        math.sin(2.0 * math.pi * (index + 0.5) / 37.0) for index in range(sample_count)
    ]
    q3 = [
        math.cos(2.0 * math.pi * (index + 0.5) / 53.0) for index in range(sample_count)
    ]
    sex = [1 if index % 2 == 0 else 2 for index in range(sample_count)]

    qt_score = polygenic_score(dosages, start=0, step=31, count=64)
    bt_score = polygenic_score(dosages, start=13, step=37, count=64)
    age_scaled = standardise([float(value) for value in age])
    qt_noise = normal_draws(rng, sample_count)
    bt_noise = normal_draws(rng, sample_count)

    quantitative = [
        0.65 * qt_score[index]
        + 0.18 * age_scaled[index]
        + 0.12 * q2[index]
        + 0.55 * qt_noise[index]
        for index in range(sample_count)
    ]
    liability = [
        0.60 * bt_score[index]
        + 0.15 * q1[index]
        - 0.12 * q3[index]
        + 0.60 * bt_noise[index]
        for index in range(sample_count)
    ]
    case_indices = {
        index
        for index, _ in sorted(
            enumerate(liability), key=lambda item: (-item[1], item[0])
        )[:cases]
    }
    binary = [2 if index in case_indices else 1 for index in range(sample_count)]

    write_tsv(
        output_dir / f"{prefix}.pheno",
        ["FID", "IID", "QT", "BT"],
        [
            [sample, sample, f"{quantitative[index]:.6f}", binary[index]]
            for index, sample in enumerate(samples)
        ],
    )
    write_tsv(
        output_dir / f"{prefix}.qcovar",
        ["FID", "IID", "AGE", "Q1", "Q2", "Q3"],
        [
            [
                sample,
                sample,
                age[index],
                f"{q1[index]:.6f}",
                f"{q2[index]:.6f}",
                f"{q3[index]:.6f}",
            ]
            for index, sample in enumerate(samples)
        ],
    )
    write_tsv(
        output_dir / f"{prefix}.catcovar",
        ["FID", "IID", "SEX"],
        [[sample, sample, sex[index]] for index, sample in enumerate(samples)],
    )


def main() -> None:
    args = parse_args()
    chromosomes = [int(value) for value in args.chromosomes.split(",")]
    if chromosomes != sorted(set(chromosomes)):
        raise ValueError("chromosomes must be unique and sorted")
    if any(chromosome < 1 or chromosome > 22 for chromosome in chromosomes):
        raise ValueError("only autosomes 1 through 22 are supported")
    if args.samples < 2:
        raise ValueError("at least two samples are required")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    samples = [f"S{index:03d}" for index in range(1, args.samples + 1)]
    rng = random.Random(args.seed)
    records, dosages = build_genotypes(
        rng, chromosomes, samples, args.variants_per_chromosome
    )
    write_vcf(args.output_dir / f"{args.prefix}_all.vcf", samples, records)
    write_sidecars(args.output_dir, args.prefix, samples, dosages, args.cases, rng)


if __name__ == "__main__":
    main()
