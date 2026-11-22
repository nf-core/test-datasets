# Files were downloaded using the following links on the 29th of August 2024

## mirgenedb_mature  = "https://mirgenedb.org/fasta/hsa?mat=1"
## mirgenedb_hairpin = "https://mirgenedb.org/static/data/hsa/hsa-pre.fas"
## mirgenedb_gff     = "https://mirgenedb.org/gff/hsa?sort=pos&all=1"

### Files were ingested instead of using them directly from the link, because they get downloaded with unusual names, and nextflow can't stage them properly

## mirgenedb_star    = derived from mirgenedb_hsa_hairpin.fa and mirgenedb_hsa.gff (see below)

`mirgenedb_hsa_star.fa` contains the 503 human star-arm sequences (miRNAs whose GFF `ID` ends with `*`, e.g. `Hsa-Mir-3613_3p*`) of the MirGeneDB v2.0 files above. It is used by the `--mirgenedb_star` test of nf-core/smrnaseq.

The live "Star sequences" download (`https://mirgenedb.org/fasta/hsa?star=1`) is now MirGeneDB v3.0 and does not match the v2.0 GFF and precursor files above. The star sequences were therefore cut out of `mirgenedb_hsa_hairpin.fa` at the GFF coordinates of each star miRNA, relative to the precursor (`pre_miRNA`) that contains it. The same derivation applied to the non-star miRNAs reproduces all 630 records of `mirgenedb_hsa_mature.fa` exactly.

```bash
python3 derive_arms.py mirgenedb_hsa.gff mirgenedb_hsa_hairpin.fa star > mirgenedb_hsa_star.fa
# check: python3 derive_arms.py mirgenedb_hsa.gff mirgenedb_hsa_hairpin.fa mature  ==  mirgenedb_hsa_mature.fa
```

`derive_arms.py`:

```python
#!/usr/bin/env python3
"""Cut MirGeneDB miRNA arm sequences out of the precursor FASTA at their GFF coordinates.

Usage: python3 derive_arms.py <mirgenedb.gff> <hairpin.fa> star|mature > out.fa
  star   -> miRNAs whose ID ends with '*' (star arms)
  mature -> all other miRNAs (dominant arms)
"""
import sys

gff_path, hairpin_path, which = sys.argv[1], sys.argv[2], sys.argv[3]

pre_seq, name = {}, None
for line in open(hairpin_path):
    line = line.strip()
    if line.startswith(">"):
        name = line[1:].split()[0]
        pre_seq[name] = ""
    elif line:
        pre_seq[name] += line

pres, arms = [], []
for line in open(gff_path):
    if line.startswith("#") or not line.strip():
        continue
    f = line.rstrip("\n").split("\t")
    attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
    rec = (f[0], int(f[3]), int(f[4]), f[6], attrs["ID"])
    if f[2] == "pre_miRNA":
        pres.append(rec)
    elif f[2] == "miRNA" and attrs["ID"].endswith("*") == (which == "star"):
        arms.append(rec)

for chrom, start, end, strand, arm in arms:
    # precursor on the same chromosome and strand whose coordinates contain the arm
    pre = next(p for p in pres if p[0] == chrom and p[3] == strand and p[1] <= start and p[2] >= end and p[4] in pre_seq)
    offset = start - pre[1] if strand == "+" else pre[2] - end
    print(f">{arm}\n{pre_seq[pre[4]][offset:offset + end - start + 1]}")
```
