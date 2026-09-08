#!/usr/bin/env python3
"""
Rebuild the CADD prescored test fixtures for the nf-core/raredisease ``annotate_cadd``
subworkflow against the minimal 9-region GIAB dataset (``raredisease`` branch).

Outputs (written under ``--repo``, default = this script's parent repo):

  resources_remapped/cadd/cadd_prescored.tar.gz
      prescored/GRCh38_v<version>/no_anno/grch38_cadd_snvs.tsv.gz     (+ .tbi)
      prescored/GRCh38_v<version>/no_anno/grch38_cadd_indels.tsv.gz   (+ .tbi)
  resources_remapped/cadd/cadd_annotations.tar.gz
      annotations/.gitkeep
  subworkflow_fixtures/cadd_test_variants.vcf.gz                      (+ .tbi)

Why this exists
---------------
``annotate_cadd`` renames the reference contigs to CADD's no-``chr`` convention
before scoring and back afterwards, and the nf-test only exercises the
all-prescored path: every test variant is present in the prescored files, so
``CADD.sh`` never touches the ~200 GB annotation bundle or a reference genome.
The fixture therefore only needs a handful of real CADD rows, lifted onto the
sliced local coordinates via ``manifests/coordinate_provenance.tsv``.

How to run
----------
Runs from any working directory (paths are resolved relative to the script, not
the cwd) with any Python >= 3.8:

  python3 scripts/build_cadd_prescored_fixture.py --snv-prescored WGS_SNVS.tsv.gz --indel-prescored GNOMAD_INDELS.tsv.gz

Takes a few seconds. Writes straight into the working tree; review with
``git diff`` / ``git status`` and commit.

Inputs
------
  * ``tabix`` and ``bgzip`` (htslib) must be on PATH.
  * Read automatically from this repo (via ``--repo``, default = the repo the
    script lives in):
      manifests/coordinate_provenance.tsv
      reference_sliced/minimal_reference.fasta
  * Passed on the command line -- the full CADD prescored sets for the target
    version. Too big to keep in the repo; download once and keep locally. The
    script only reads them locally through ``tabix`` (no network access), so
    each needs its ``.tbi`` next to it.
      --snv-prescored    whole_genome_SNVs.tsv.gz       ~80 GB, chr-prefixed
      --indel-prescored  gnomad.genomes.r4.0.indel.tsv.gz ~1 GB, no chr prefix
    Both from the CADD release tree (pick the folder for the target version):
      https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz          (+ .tbi)
      https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz  (+ .tbi)
    Landing page / other versions: https://cadd.gs.washington.edu/download

Reproduce the current (v1.7) fixture
-----------------------------------
  python3 scripts/build_cadd_prescored_fixture.py --snv-prescored whole_genome_SNVs.tsv.gz --indel-prescored gnomad.genomes.r4.0.indel.tsv.gz

Next CADD update
----------------
  python3 scripts/build_cadd_prescored_fixture.py --cadd-version 1.8 --snv-prescored WGS_SNVS.tsv.gz --indel-prescored GNOMAD_INDELS.tsv.gz

Then, in nf-core/raredisease:
  * bump the prescored path version in
    subworkflows/local/annotate_cadd/tests/nextflow_real.config only if CADD.sh's
    ``config_GRCh38_v<version>_noanno.yml`` -> ``PrescoredFolder`` changed, and
  * re-record subworkflows/local/annotate_cadd/tests/main.nf.test.snap
    (nf-test ... --update-snapshot).
"""

from __future__ import annotations

import argparse
import random
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path
from typing import NamedTuple

# Nuclear sliced contigs, in coordinate_provenance.tsv order. chrM is skipped:
# CADD publishes no mitochondrial prescored scores.
NUCLEAR_CONTIGS = ["chr1", "chr7", "chr12", "chr16", "chr20", "chr21", "chrX"]

# Column indices in a CADD prescored TSV row: Chrom Pos Ref Alt RawScore PHRED
CADD_COL_POS, CADD_COL_REF, CADD_COL_ALT, CADD_COL_RAW, CADD_COL_PHRED = 1, 2, 3, 4, 5

END_MARGIN_BP = 200          # keep picks at least this far from either contig end
PICKS_PER_CONTIG = 3         # SNVs and indels, each, per contig
DEFAULT_SEED = 42

CADD_LICENCE_LINE = (
    "## CADD GRCh38-v{version} (c) University of Washington, Hudson-Alpha Institute "
    "for Biotechnology and Berlin Institute of Health at Charite - "
    "Universitaetsmedizin Berlin 2013-2023. All rights reserved."
)
CADD_COLUMN_HEADER = "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED"


class ProvenanceSegment(NamedTuple):
    """One contiguous real->local coordinate mapping from coordinate_provenance.tsv."""

    local_start: int
    local_end: int
    real_chrom: str
    real_start: int
    real_end: int

    @property
    def real_to_local_offset(self) -> int:
        return self.local_start - self.real_start


class PrescoredVariant(NamedTuple):
    """A CADD prescored row already lifted onto sliced local coordinates."""

    contig: str          # no "chr" prefix, e.g. "20"
    pos: int             # 1-based, sliced local coordinate
    ref: str
    alt: str
    raw_score: str
    phred: str


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def run_checked(command: list[str], **kwargs) -> subprocess.CompletedProcess:
    return subprocess.run(command, check=True, **kwargs)


def strip_chr(contig: str) -> str:
    return contig[3:] if contig.startswith("chr") else contig


def contig_sort_key(contig_without_chr: str) -> int:
    rank = {str(number): number for number in range(1, 23)}
    rank.update({"X": 23, "Y": 24, "M": 25})
    return rank.get(contig_without_chr, 99)


# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #
def load_provenance(provenance_tsv: Path) -> dict[str, list[ProvenanceSegment]]:
    """Map each nuclear contig to its provenance segments (chr20 has two, the
    rest have one)."""
    segments_by_contig: dict[str, list[ProvenanceSegment]] = {}
    with provenance_tsv.open() as tsv:
        header_fields = tsv.readline().rstrip("\n").split("\t")
        column_index = {name: index for index, name in enumerate(header_fields)}
        for line in tsv:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            local_contig = fields[column_index["new_contig"]]
            if local_contig not in NUCLEAR_CONTIGS:
                continue
            segments_by_contig.setdefault(local_contig, []).append(
                ProvenanceSegment(
                    local_start=int(fields[column_index["new_start"]]),
                    local_end=int(fields[column_index["new_end"]]),
                    real_chrom=fields[column_index["original_chrom"]],
                    real_start=int(fields[column_index["original_start"]]),
                    real_end=int(fields[column_index["original_end"]]),
                )
            )
    missing = [contig for contig in NUCLEAR_CONTIGS if contig not in segments_by_contig]
    if missing:
        raise SystemExit(f"coordinate_provenance.tsv has no rows for: {', '.join(missing)}")
    return segments_by_contig


def load_reference(reference_fasta: Path) -> dict[str, str]:
    """Contig (as written in the FASTA, e.g. "chr1") -> uppercase sequence.
    Index with ``pos - 1`` for 1-based coordinates."""
    sequence_chunks: dict[str, list[str]] = {}
    current_contig = None
    with reference_fasta.open() as fasta:
        for line in fasta:
            if line.startswith(">"):
                current_contig = line[1:].split()[0]
                sequence_chunks[current_contig] = []
            else:
                sequence_chunks[current_contig].append(line.strip())
    return {contig: "".join(chunks).upper() for contig, chunks in sequence_chunks.items()}


# --------------------------------------------------------------------------- #
# Extract + remap
# --------------------------------------------------------------------------- #
def remap_prescored_contig(
    prescored_gz: Path,
    local_contig: str,
    segments: list[ProvenanceSegment],
    file_is_chr_prefixed: bool,
):
    """tabix every provenance segment of ``local_contig`` out of ``prescored_gz``
    and yield each row as a :class:`PrescoredVariant` on sliced local
    coordinates. Rows that fall outside the mapped local window are dropped."""
    local_contig_nochr = strip_chr(local_contig)
    for segment in segments:
        query_chrom = segment.real_chrom if file_is_chr_prefixed else strip_chr(segment.real_chrom)
        region = f"{query_chrom}:{segment.real_start}-{segment.real_end}"
        tabix_result = run_checked(
            ["tabix", str(prescored_gz), region],
            stdout=subprocess.PIPE,
            text=True,
        )
        for line in tabix_result.stdout.splitlines():
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            local_pos = int(fields[CADD_COL_POS]) + segment.real_to_local_offset
            if not segment.local_start <= local_pos <= segment.local_end:
                continue
            yield PrescoredVariant(
                contig=local_contig_nochr,
                pos=local_pos,
                ref=fields[CADD_COL_REF],
                alt=fields[CADD_COL_ALT],
                raw_score=fields[CADD_COL_RAW],
                phred=fields[CADD_COL_PHRED],
            )


def collect_prescored_variants(
    prescored_gz: Path,
    provenance: dict[str, list[ProvenanceSegment]],
    file_is_chr_prefixed: bool,
) -> list[PrescoredVariant]:
    variants: list[PrescoredVariant] = []
    for local_contig in NUCLEAR_CONTIGS:
        count_before = len(variants)
        variants.extend(
            remap_prescored_contig(
                prescored_gz, local_contig, provenance[local_contig], file_is_chr_prefixed
            )
        )
        log(f"  {local_contig}: {len(variants) - count_before} prescored rows")
    return variants


# --------------------------------------------------------------------------- #
# Pick the fixture variants
# --------------------------------------------------------------------------- #
def ref_matches_reference(reference: dict[str, str], contig_nochr: str, pos: int, ref: str) -> bool:
    sequence = reference.get("chr" + contig_nochr) or reference.get(contig_nochr)
    if sequence is None or pos < 1 or pos - 1 + len(ref) > len(sequence):
        return False
    return sequence[pos - 1 : pos - 1 + len(ref)] == ref.upper()


def pick_fixture_variants(
    candidates: list[PrescoredVariant],
    reference: dict[str, str],
    random_generator: random.Random,
    kind: str,
) -> list[PrescoredVariant]:
    """Choose ``PICKS_PER_CONTIG`` distinct positions per contig, at least
    ``END_MARGIN_BP`` from either end, with REF verified against the sliced
    reference. Deterministic under the seed: every population handed to the RNG
    is a sorted list."""
    picked: list[PrescoredVariant] = []
    for local_contig in NUCLEAR_CONTIGS:
        contig_nochr = strip_chr(local_contig)
        contig_length = len(reference.get(local_contig) or reference.get(contig_nochr))
        candidates_by_pos: dict[int, list[PrescoredVariant]] = {}
        for variant in candidates:
            if variant.contig != contig_nochr:
                continue
            if variant.pos <= END_MARGIN_BP or variant.pos > contig_length - END_MARGIN_BP:
                continue
            if not ref_matches_reference(reference, contig_nochr, variant.pos, variant.ref):
                continue
            candidates_by_pos.setdefault(variant.pos, []).append(variant)

        usable_positions = sorted(candidates_by_pos)
        if len(usable_positions) < PICKS_PER_CONTIG:
            raise SystemExit(
                f"{kind}: only {len(usable_positions)} usable positions on {local_contig}"
            )
        for pos in sorted(random_generator.sample(usable_positions, PICKS_PER_CONTIG)):
            alleles_at_pos = sorted(candidates_by_pos[pos], key=lambda candidate: (candidate.ref, candidate.alt))
            picked.append(random_generator.choice(alleles_at_pos))
    return picked


# --------------------------------------------------------------------------- #
# Write outputs
# --------------------------------------------------------------------------- #
def write_prescored_tsv(output_gz: Path, variants: list[PrescoredVariant], cadd_version: str) -> None:
    output_gz.parent.mkdir(parents=True, exist_ok=True)
    uncompressed = output_gz.with_suffix("")  # drop the .gz
    ordered = sorted(variants, key=lambda variant: (contig_sort_key(variant.contig), variant.pos, variant.ref, variant.alt))
    with uncompressed.open("w") as tsv:
        tsv.write(CADD_LICENCE_LINE.format(version=cadd_version) + "\n")
        tsv.write(CADD_COLUMN_HEADER + "\n")
        for variant in ordered:
            tsv.write(
                "\t".join(
                    [variant.contig, str(variant.pos), variant.ref, variant.alt,
                     variant.raw_score, variant.phred]
                )
                + "\n"
            )
    run_checked(["bgzip", "-f", str(uncompressed)])
    run_checked(["tabix", "-f", "-s", "1", "-b", "2", "-e", "2", "-c", "#", str(output_gz)])


def write_test_vcf(
    output_gz: Path,
    snvs: list[PrescoredVariant],
    indels: list[PrescoredVariant],
    reference: dict[str, str],
) -> None:
    """chr-prefixed local-coordinate VCF with the same 42 sites, for the
    subworkflow's ``ch_vcf`` input (annotate_cadd strips the prefix itself)."""
    output_gz.parent.mkdir(parents=True, exist_ok=True)
    uncompressed = output_gz.with_suffix("")
    all_variants = sorted(
        [*snvs, *indels], key=lambda variant: (contig_sort_key(variant.contig), variant.pos, variant.ref, variant.alt)
    )
    with uncompressed.open("w") as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        for local_contig in NUCLEAR_CONTIGS:
            contig_length = len(reference.get(local_contig) or reference.get(strip_chr(local_contig)))
            vcf.write(f"##contig=<ID={local_contig},length={contig_length}>\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for variant in all_variants:
            vcf.write(f"chr{variant.contig}\t{variant.pos}\t.\t{variant.ref}\t{variant.alt}\t.\t.\t.\n")
    run_checked(["bgzip", "-f", str(uncompressed)])
    run_checked(["tabix", "-f", "-p", "vcf", str(output_gz)])


def make_tarball(output_tarball: Path, staging_dir: Path, top_level_dir: str) -> None:
    with tarfile.open(output_tarball, "w:gz") as tarball:
        tarball.add(staging_dir / top_level_dir, arcname=top_level_dir)


# --------------------------------------------------------------------------- #
def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    default_repo = Path(__file__).resolve().parents[1]
    parser.add_argument(
        "--repo", type=Path, default=default_repo,
        help=f"test-datasets checkout (default: {default_repo})",
    )
    parser.add_argument(
        "--snv-prescored", type=Path, required=True,
        help="CADD whole_genome_SNVs.tsv.gz (bgzipped + .tbi)",
    )
    parser.add_argument(
        "--indel-prescored", type=Path, required=True,
        help="CADD gnomAD indel tsv.gz (bgzipped + .tbi)",
    )
    parser.add_argument(
        "--cadd-version", default="1.7",
        help="CADD version string used in the prescored path and licence line (default: 1.7)",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()

    for tool in ("tabix", "bgzip"):
        if shutil.which(tool) is None:
            raise SystemExit(f"{tool} not found on PATH")

    provenance_tsv = args.repo / "manifests" / "coordinate_provenance.tsv"
    reference_fasta = args.repo / "reference_sliced" / "minimal_reference.fasta"
    for required_input in (provenance_tsv, reference_fasta, args.snv_prescored, args.indel_prescored):
        if not required_input.exists():
            raise SystemExit(f"missing input: {required_input}")

    random_generator = random.Random(args.seed)

    log("Loading provenance + sliced reference ...")
    provenance = load_provenance(provenance_tsv)
    reference = load_reference(reference_fasta)

    log("Extracting + remapping SNV prescored rows ...")
    snv_candidates = collect_prescored_variants(args.snv_prescored, provenance, file_is_chr_prefixed=True)
    log("Extracting + remapping indel prescored rows ...")
    indel_candidates = collect_prescored_variants(args.indel_prescored, provenance, file_is_chr_prefixed=False)

    log("Picking fixture variants ...")
    snv_picks = pick_fixture_variants(snv_candidates, reference, random_generator, "SNV")
    indel_picks = pick_fixture_variants(indel_candidates, reference, random_generator, "indel")
    log(f"  {len(snv_picks)} SNVs + {len(indel_picks)} indels")

    cadd_dir = args.repo / "resources_remapped" / "cadd"
    cadd_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp_name:
        staging_dir = Path(tmp_name)
        prescored_dir = staging_dir / "prescored" / f"GRCh38_v{args.cadd_version}" / "no_anno"
        write_prescored_tsv(prescored_dir / "grch38_cadd_snvs.tsv.gz", snv_picks, args.cadd_version)
        write_prescored_tsv(prescored_dir / "grch38_cadd_indels.tsv.gz", indel_picks, args.cadd_version)
        (staging_dir / "annotations").mkdir()
        (staging_dir / "annotations" / ".gitkeep").touch()

        make_tarball(cadd_dir / "cadd_prescored.tar.gz", staging_dir, "prescored")
        make_tarball(cadd_dir / "cadd_annotations.tar.gz", staging_dir, "annotations")

    write_test_vcf(
        args.repo / "subworkflow_fixtures" / "cadd_test_variants.vcf.gz",
        snv_picks,
        indel_picks,
        reference,
    )

    log("\nWrote:")
    for relative_path in (
        "resources_remapped/cadd/cadd_prescored.tar.gz",
        "resources_remapped/cadd/cadd_annotations.tar.gz",
        "subworkflow_fixtures/cadd_test_variants.vcf.gz",
        "subworkflow_fixtures/cadd_test_variants.vcf.gz.tbi",
    ):
        log(f"  {relative_path}")


if __name__ == "__main__":
    main()
