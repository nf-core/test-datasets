"""Cut the chrX test reference down to the regions around a few expressed genes.

Usage: make_chrx_subset.py X.fa genes_chrX.gtf read_alignments.tsv n_genes keep_genes n_pairs out_dir

read_alignments.tsv holds `samtools view -F 0x904` columns 1-8 of every sample
(name, flag, chrom, pos, mapq, cigar, rnext, pnext).
n_genes genes with at least two transcripts are taken by their number of reads, and
keep_genes, a comma separated list of gene IDs, is kept whatever its counts.
n_pairs is the number of read pairs kept per sample, drawn at random with a fixed seed.

Writes out_dir/X.fa, out_dir/genes_chrX.gtf, out_dir/regions.tsv and
out_dir/keep_reads.txt, the names of the pairs to keep. The regions are joined into
one contig named X, separated by runs of N, so the contig name stays the same and
only the coordinates change. The GTF keeps the features and attributes the pipeline
uses, to keep it small.
"""

import bisect
import collections
import random
import re
import sys

fasta, gtf, aln, n_genes, keep_genes, n_pairs, out_dir = sys.argv[1:]
n_genes, n_pairs = int(n_genes), int(n_pairs)
keep_genes = [g for g in keep_genes.split(",") if g]
PAD = 1000
GAP = 10000
SPACER = 1000
FEATURES = {"transcript", "exon", "CDS"}
ATTRIBUTES = ["gene_id", "transcript_id", "gene_name", "gene_biotype", "exon_number"]

# Gene spans and transcripts from the exons
genes = {}
transcripts = collections.defaultdict(set)
for line in open(gtf):
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 9 or fields[2] != "exon":
        continue
    gene = re.search(r'gene_id "([^"]+)"', fields[8]).group(1)
    transcripts[gene].add(re.search(r'transcript_id "([^"]+)"', fields[8]).group(1))
    start, end = int(fields[3]), int(fields[4])
    s, e = genes.get(gene, (start, end))
    genes[gene] = (min(s, start), max(e, end))

# Reads per gene, from the alignment starts of all samples
starts = []
mates = {}
for line in open(aln):
    name, flag, chrom, pos, _mapq, _cigar, rnext, pnext = line.rstrip("\n").split("\t")[:8]
    starts.append(int(pos))
    if int(flag) & 64 and rnext == "=":
        mates[name] = (int(pos), int(pnext))
starts.sort()
counts = {g: bisect.bisect_right(starts, e) - bisect.bisect_left(starts, s) for g, (s, e) in genes.items()}
ranked = [g for g in sorted(counts, key=lambda g: -counts[g]) if len(transcripts[g]) >= 2]
selected = set(ranked[:n_genes]) | set(keep_genes)

# Pad and merge the selected genes, then widen each region until no gene crosses its ends
spans = sorted(genes.values())
regions = sorted([genes[g][0] - PAD, genes[g][1] + PAD] for g in selected)
changed = True
while changed:
    changed = False
    merged = []
    for s, e in regions:
        if merged and s <= merged[-1][1] + GAP:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    for region in merged:
        for gs, ge in spans:
            if gs <= region[1] and ge >= region[0] and (gs < region[0] or ge > region[1]):
                region[0], region[1] = min(region[0], gs), max(region[1], ge)
                changed = True
    regions = merged

seq = "".join(line.strip() for line in open(fasta) if not line.startswith(">"))
regions = [[max(1, s), min(len(seq), e)] for s, e in regions]

# New contig and coordinate offsets
parts, offsets, position = [], [], 1
for s, e in regions:
    offsets.append((s, e, position - s))
    parts.append(seq[s - 1:e])
    position += e - s + 1
    parts.append("N" * SPACER)
    position += SPACER
contig = "".join(parts[:-1])

with open(f"{out_dir}/X.fa", "w") as handle:
    handle.write(">X\n")
    for i in range(0, len(contig), 60):
        handle.write(contig[i:i + 60] + "\n")

with open(f"{out_dir}/regions.tsv", "w") as handle:
    handle.write("old_start\told_end\tnew_start\tnew_end\n")
    for s, e, shift in offsets:
        handle.write(f"{s}\t{e}\t{s + shift}\t{e + shift}\n")

region_starts = [s for s, _, _ in offsets]


def shift_of(start, end):
    i = bisect.bisect_right(region_starts, start) - 1
    if i >= 0 and offsets[i][1] >= end:
        return offsets[i][2]
    return None


# Annotation of the genes inside the regions, in the new coordinates
kept_genes = {g for g, (s, e) in genes.items() if shift_of(s, e) is not None}
with open(f"{out_dir}/genes_chrX.gtf", "w") as handle:
    for line in open(gtf):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] not in FEATURES:
            continue
        attributes = dict(re.findall(r'(\S+) "([^"]*)"', fields[8]))
        if attributes["gene_id"] not in kept_genes:
            continue
        shift = shift_of(int(fields[3]), int(fields[4]))
        fields[3], fields[4] = str(int(fields[3]) + shift), str(int(fields[4]) + shift)
        fields[8] = " ".join(f'{key} "{attributes[key]}";' for key in ATTRIBUTES if key in attributes)
        handle.write("\t".join(fields) + "\n")

# Pairs whose two mates start within the same region, at most n_pairs per sample
by_sample = collections.defaultdict(list)
for name, (pos, mate_pos) in mates.items():
    a, b = shift_of(pos, pos), shift_of(mate_pos, mate_pos)
    if a is not None and a == b:
        by_sample[name.split(".")[0]].append(name)
random.seed(100)
with open(f"{out_dir}/keep_reads.txt", "w") as handle:
    for sample in sorted(by_sample):
        names = sorted(by_sample[sample])
        for name in sorted(random.sample(names, min(n_pairs, len(names)))):
            handle.write(name + "\n")

print(
    f"{len(regions)} regions, {len(contig) / 1e6:.2f} Mb, {len(kept_genes)} genes (of {len(genes)}), "
    + ", ".join(f"{s} {min(n_pairs, len(n))} of {len(n)} pairs" for s, n in sorted(by_sample.items())),
    file=sys.stderr,
)
