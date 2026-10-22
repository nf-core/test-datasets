"""Turn an unstranded paired-end sample into a reverse stranded one.

Usage: make_reverse_stranded.py genes.gtf read1_alignments.tsv in_1.fq.gz in_2.fq.gz out_1.fq.gz out_2.fq.gz

read1_alignments.tsv holds `samtools view -f 0x40 -F 0x904` columns 1-6 (name, flag, chrom, pos, mapq, cigar).
A pair whose read 1 lies on the strand of the transcripts it overlaps comes from the sense strand;
swapping its mates makes read 1 antisense, as in a dUTP (reverse stranded) library.
Pairs that overlap transcripts on both strands, or no transcript, are left as they are.
"""

import bisect
import gzip
import re
import sys

gtf, aln, in1, in2, out1, out2 = sys.argv[1:]

# Transcript spans per chromosome and strand
transcripts = {}
with open(gtf) as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if fields[2] != "transcript":
            continue
        transcripts.setdefault((fields[0], fields[6]), []).append((int(fields[3]), int(fields[4])))

index = {}
for key, spans in transcripts.items():
    spans.sort()
    starts = [s for s, _ in spans]
    # Running maximum end, so a binary search on starts finds every overlapping span
    max_ends, current = [], 0
    for _, end in spans:
        current = max(current, end)
        max_ends.append(current)
    index[key] = (starts, spans, max_ends)


def overlaps(chrom, strand, start, end):
    if (chrom, strand) not in index:
        return False
    starts, spans, max_ends = index[(chrom, strand)]
    i = bisect.bisect_right(starts, end) - 1
    while i >= 0 and max_ends[i] >= start:
        if spans[i][1] >= start:
            return True
        i -= 1
    return False


swap = set()
counts = {"sense": 0, "antisense": 0, "ambiguous": 0, "intergenic": 0}
with open(aln) as handle:
    for line in handle:
        name, flag, chrom, pos, _mapq, cigar = line.rstrip("\n").split("\t")[:6]
        start = int(pos)
        end = start + sum(int(n) for n, op in re.findall(r"(\d+)([MDN=X])", cigar)) - 1
        read_strand = "-" if int(flag) & 16 else "+"
        on_plus = overlaps(chrom, "+", start, end)
        on_minus = overlaps(chrom, "-", start, end)
        if on_plus and on_minus:
            counts["ambiguous"] += 1
        elif not (on_plus or on_minus):
            counts["intergenic"] += 1
        elif (read_strand == "+") == on_plus:
            counts["sense"] += 1
            swap.add(name)
        else:
            counts["antisense"] += 1


def records(path):
    with gzip.open(path, "rt") as handle:
        while True:
            record = [handle.readline() for _ in range(4)]
            if not record[0]:
                return
            yield record


swapped = 0
with gzip.open(out1, "wt", compresslevel=9) as o1, gzip.open(out2, "wt", compresslevel=9) as o2:
    for r1, r2 in zip(records(in1), records(in2)):
        name = r1[0][1:].split()[0]
        assert name == r2[0][1:].split()[0], f"Mates out of order at {name}"
        if name in swap:
            # Keep the /1 and /2 headers, swap sequences and qualities
            r1, r2 = [r1[0], r2[1], r1[2], r2[3]], [r2[0], r1[1], r2[2], r1[3]]
            swapped += 1
        o1.writelines(r1)
        o2.writelines(r2)

print(out1, counts, "swapped", swapped, file=sys.stderr)
