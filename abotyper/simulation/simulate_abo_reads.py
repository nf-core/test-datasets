#!/usr/bin/env python3
"""Simulate nanopore ABO amplicon reads for the nf-core/abotyper test data.

No human reads are used. Each allele sequence is the exon 2-7 amplicon of the ABO RefSeqGene NG_006669.2 with the
coding changes of that allele from the ISBT ABO allele table applied (introns stay as NG_006669.2). Each sample is a
pair of alleles; reads are drawn from the two alleles with a fixed random seed and given nanopore-like errors and
quality scores by the simple model below.

Requirements: Python >= 3.11 and abotyper-core (used only to place ISBT cDNA changes on NG_006669.2, with the same
exon offsets the pipeline uses).

    python simulate_abo_reads.py --ref NG_006669v2.fasta --amplicon 17922-24675 --outdir .

writes, for every sample in SAMPLES:
    <sample>.fastq.gz          simulated reads (byte-identical between runs: fixed seed, gzip mtime 0)
and, once:
    simulation/alleles.fasta   the allele amplicon sequences the reads were drawn from
    simulation/truth.tsv       sample, alleles, genotype, reads simulated per allele

Error model (per read base, before quality assignment), chosen to resemble R10.4.1 super-accuracy amplicon reads
(modal read accuracy about 99%):
    substitution 0.4%, insertion 0.2%, deletion 0.3%
    homopolymer runs of >= 4 bases: the run length changes by -1 (6%) or +1 (2%) per read
    read ends: 0-25 bases missing at each end (incomplete primer/adapter trimming is not modelled)
    strand: half of the reads are reverse-complemented
Quality scores are drawn per base around Q25 (Q7 at erroneous bases), as dorado writes them.
"""
from __future__ import annotations

import argparse
import gzip
import random
import uuid
from pathlib import Path

from abotyper_core.panel import default_panel_path, load_panel
from abotyper_core.readtyper import ISBT_TABLE, _CHANGE, exon_offsets, genomic_of, read_fasta

SEED = 20261117
# sample name (as on the abotyper test-data branch) -> (allele 1, allele 2, number of reads)
SAMPLES = {
    "ABO_AA_barcode01": ("ABO*A1.01", "ABO*A1.01", 100),
    "ABO_AB_barcode02": ("ABO*A1.01", "ABO*B.01", 200),
    "ABO_AO1_barcode03": ("ABO*A1.01", "ABO*O.01.01", 100),
    "ABO_AO2_barcode04": ("ABO*A1.01", "ABO*O.02.01", 100),
    "ABO_BB_barcode05": ("ABO*B.01", "ABO*B.01", 100),
    "ABO_BO1_barcode06": ("ABO*B.01", "ABO*O.01.01", 100),
    "ABO_BO2_barcode07": ("ABO*B.01", "ABO*O.02.01", 100),
    "ABO_O1O1_barcode08": ("ABO*O.01.01", "ABO*O.01.01", 100),
    "ABO_O1O2_barcode09": ("ABO*O.01.01", "ABO*O.02.01", 100),
    "ABO_O2O2_barcode10": ("ABO*O.02.01", "ABO*O.02.01", 100),
}

P_SUB, P_INS, P_DEL = 0.004, 0.002, 0.003
P_HP_DEL, P_HP_INS, HP_LEN = 0.06, 0.02, 4
MAX_END_LOSS = 25
COMP = str.maketrans("ACGT", "TGCA")


def allele_changes(name: str) -> list[str]:
    with open(ISBT_TABLE) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == name:
                return [c.strip() for c in f[2].split(";") if c.strip()]
    raise SystemExit(f"{name} not in {ISBT_TABLE}")


def build_allele(seq: str, changes: list[str], offsets: dict) -> str:
    """Apply ISBT cDNA changes to the 1-based genomic sequence, right to left so positions stay valid."""
    edits = []
    for ch in changes:
        m = _CHANGE.match(ch)
        if not m:
            raise SystemExit(f"cannot parse {ch}")
        g1 = genomic_of(int(m.group(1)), m.group(2), offsets)
        g2 = genomic_of(int(m.group(3)), m.group(4), offsets) if m.group(3) else g1
        if g1 is None:
            continue                                  # outside the exon 2-7 amplicon (e.g. exon 1)
        edits.append((g1, g2, m))
    s = list(seq)
    for g1, g2, m in sorted(edits, key=lambda e: e[0], reverse=True):
        if m.group(5):                                # SNV
            assert s[g1 - 1] == m.group(5), f"reference mismatch at g.{g1}"
            s[g1 - 1] = m.group(6)
        elif m.group(7) == "del":
            del s[g1 - 1:g2]
        elif m.group(7) == "dup":
            s[g2:g2] = s[g1 - 1:g2]
        else:                                         # insertion after g1
            s[g1:g1] = list(m.group(9))
    return "".join(s)


def mutate(rng: random.Random, template: str) -> tuple[str, str]:
    out, qual = [], []
    i, n = 0, len(template)
    while i < n:
        b = template[i]
        j = i
        while j < n and template[j] == b:
            j += 1
        run = j - i
        if run >= HP_LEN:                             # homopolymer length error
            r = rng.random()
            run += -1 if r < P_HP_DEL else (1 if r < P_HP_DEL + P_HP_INS else 0)
            for _ in range(run):
                out.append(b)
                qual.append(rng.randint(14, 30))
            i = j
            continue
        r = rng.random()
        if r < P_DEL:
            pass
        elif r < P_DEL + P_SUB:
            out.append(rng.choice([x for x in "ACGT" if x != b]))
            qual.append(rng.randint(3, 10))
        else:
            out.append(b)
            qual.append(rng.randint(18, 35))
            if r > 1 - P_INS:
                out.append(rng.choice("ACGT"))
                qual.append(rng.randint(3, 10))
        i += 1
    return "".join(out), "".join(chr(q + 33) for q in qual)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ref", required=True, type=Path, help="NG_006669.2 FASTA")
    ap.add_argument("--amplicon", required=True, help="1-based inclusive START-END of the amplicon on NG_006669.2")
    ap.add_argument("--outdir", default=Path("."), type=Path)
    a = ap.parse_args()
    start, end = (int(x) for x in a.amplicon.split("-"))

    seq = read_fasta(a.ref)
    offsets = exon_offsets(load_panel(default_panel_path()))
    names = sorted({x for s in SAMPLES.values() for x in s[:2]})
    alleles = {}
    for nm in names:
        full = build_allele(seq, allele_changes(nm), offsets)
        shift = len(full) - len(seq)                  # all changes lie inside the amplicon
        alleles[nm] = full[start - 1:end + shift]

    sim = a.outdir / "simulation"
    sim.mkdir(parents=True, exist_ok=True)
    with open(sim / "alleles.fasta", "w", newline="\n") as fh:
        for nm in names:
            s = alleles[nm]
            fh.write(f">{nm} NG_006669.2:{start}-{end} with ISBT changes {';'.join(allele_changes(nm)) or 'none'}\n")
            fh.writelines(s[k:k + 80] + "\n" for k in range(0, len(s), 80))

    truth = ["sample\tallele_1\tallele_2\tgenotype\treads_allele_1\treads_allele_2"]
    rng = random.Random(SEED)
    for sample, (a1, a2, nreads) in SAMPLES.items():
        n1 = sum(rng.random() < 0.5 for _ in range(nreads))
        picks = [a1] * n1 + [a2] * (nreads - n1)
        rng.shuffle(picks)
        records = []
        for al in picks:
            t = alleles[al]
            t = t[rng.randint(0, MAX_END_LOSS):len(t) - rng.randint(0, MAX_END_LOSS)]
            if rng.random() < 0.5:
                t = t.translate(COMP)[::-1]
            s, q = mutate(rng, t)
            rid = uuid.UUID(int=rng.getrandbits(128), version=4)
            records.append(f"@{rid} simulated\n{s}\n+\n{q}\n")
        with open(a.outdir / f"{sample}.fastq.gz", "wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
                gz.write("".join(records).encode())
        g = "".join(sorted(_class(x) for x in (a1, a2)))
        truth.append(f"{sample}\t{a1}\t{a2}\t{g}\t{n1}\t{nreads - n1}")
    with open(sim / "truth.tsv", "w", newline="\n") as fh:
        fh.write("\n".join(truth) + "\n")


def _class(name: str) -> str:
    base = name.split("*")[1]
    return "A" if base.startswith("A") else ("B" if base.startswith("B") else "O")


if __name__ == "__main__":
    main()
