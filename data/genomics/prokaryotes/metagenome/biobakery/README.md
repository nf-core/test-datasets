# bioBakery test data

Test datasets for [bioBakery](https://huttenhower.sph.harvard.edu/biobakery-workflows/) tools.

## MetaPhlAn toy database

A minimal CHOCOPhlAn database for running MetaPhlAn without downloading the full 60 GB index.
Built from `mpa_vJan25_CHOCOPhlAnSGB_202503` (produced 2026-04-30) by randomly sampling 10 SGBs
per kingdom (Bacteria, Archaea, Eukaryota) and keeping up to 200 markers per SGB, then building
a Bowtie2 large index from the resulting FASTA.

| Metric              | Value                        |
| ------------------- | ---------------------------- |
| Kingdoms            | Bacteria, Archaea, Eukaryota |
| SGBs per kingdom    | 10 (30 total)                |
| Max markers per SGB | 200 (6,000 total)            |
| Archive size        | ~14 MB                       |

Files:

- `mpa_vJan25_TOY_CHOCOPhlAnSGB_202503.tar.gz` — archive containing the marker database pickle,
  marker sequences (FASTA), and Bowtie2 large index files.
- `test_genefamilies.tsv` – mock file for gene family abundance data, to test the HUMAnN3 module.

## HUMAnN 3 demo databases

Small databases, reads and a taxonomic profile for running HUMAnN 3.6 in module tests, in
`humann3/`. With them HUMAnN finds 11 gene families (6 from the nucleotide search, 5 from the
translated search) and two MetaCyc pathways (`PWY-5423`, `PWY490-3`).

`uniref_nfDEMO.tar.gz` is copied unchanged from `data/database/humann/v3/` on the `funcprofiler`
branch, which built it from the demo databases bundled with HUMAnN 3.6.1 (`humann/data/`, folder
renamed to `uniref_nfDEMO` so HUMAnN does not warn about demo databases; commands in the
`funcprofiler` branch `README.md`). `chocophlan_nfDEMO.tar.gz` and `utility_nfDEMO.tar.gz` are the
`funcprofiler` ones cut further, as below.

`demo_profile.txt` is a hand-written MetaPhlAn 4 profile with a `vJan21` header (HUMAnN 3.6 accepts
only MetaPhlAn 3 and `vJan21` profiles), listing _Bacteroides dorei_ and _Bacteroides vulgatus_ at 50% each, the two
species of the demo nucleotide database.

`demo_subset.fastq.gz` (177 reads) and `chocophlan_nfDEMO.tar.gz` (three genes per species) come from
one HUMAnN 3.6.1 run on the HUMAnN demo reads (`examples/demo.fastq.gz`, 21,000 reads) with the full
`funcprofiler` nucleotide and protein databases:

```bash
# Inputs: chocophlan/ and uniref/ = the extracted funcprofiler chocophlan_nfDEMO.tar.gz and
# uniref_nfDEMO.tar.gz; demo.fastq.gz = examples/demo.fastq.gz of HUMAnN 3.6.1; demo_profile.txt.
# Run with HUMAnN 3.6.1 on the PATH.
humann --threads 4 --input demo.fastq.gz --taxonomic-profile demo_profile.txt \
    --nucleotide-database chocophlan --protein-database uniref --output-basename demo --output full
T=full/demo_humann_temp

# Genes: per species, the three gene families with the highest abundance, i.e. genes that pass
# the nucleotide coverage filters with their reads.
for sp in Bacteroides_dorei Bacteroides_vulgatus; do
    grep -P "\|g__Bacteroides\.s__${sp}\t" full/demo_genefamilies.tsv | sort -t$'\t' -k2,2gr | head -3 |
        cut -f1 | sed 's/|.*//' | while read -r u; do
            awk -F'\t' -v u="|${u}|" -v s="s__${sp}|" 'index($2, u) && index($2, s) { print $2 }' \
                "$T/demo_bowtie2_aligned.tsv" | sort | uniq -c | sort -k1,1nr | head -1 | awk '{ print $2 }'
        done
done > genes.txt

# Reads: those aligned to the kept genes, those DIAMOND assigns to the five proteins behind
# PWY-5423 and PWY490-3, and the unclassified reads that match nothing.
awk -F'\t' 'NR == FNR { k[$1]; next } ($2 in k) { print $1 }' genes.txt "$T/demo_bowtie2_aligned.tsv" > reads.txt
awk -F'\t' '{ split($2, a, "|"); if (a[1] ~ /^UniRef90_(Q9ZUH4|W8YTG4|V5V2L3|Z9JXD8|W7DVV5)$/) { split($1, b, "|"); print b[1] } }' \
    "$T/demo_diamond_aligned.tsv" >> reads.txt
grep '>' "$T/demo_diamond_unaligned.fa" | sed 's/^>//; s/|.*//' | grep unclassified >> reads.txt
sort -u reads.txt > reads.sorted.txt

zcat demo.fastq.gz | awk 'NR == FNR { k["@" $1]; next } FNR % 4 == 1 { p = ($1 in k) } p' reads.sorted.txt - |
    gzip -n --best > demo_subset.fastq.gz
mkdir chocophlan_min
for f in chocophlan/*.ffn.gz; do
    zcat "$f" | awk 'NR == FNR { k[">" $1]; next } /^>/ { p = ($1 in k) } p' genes.txt - |
        gzip -n --best > "chocophlan_min/$(basename "$f")"
done
tar czf chocophlan_nfDEMO.tar.gz --owner=0 --group=0 --sort=name -C chocophlan_min .
```

`utility_nfDEMO.tar.gz` is the `funcprofiler` utility database cut further: every mapping and name
file keeps only its first line, and the UniRef mapping files keep only the first three UniRef IDs of
that line. The first line of `map_level4ec_uniref90.txt.gz` (EC `1.1.1.1`) holds the three gene
families in `test_genefamilies.tsv`, so `humann_regroup_table` gives the same output as with the
larger file.

```bash
# utility/ = extracted funcprofiler utility_nfDEMO.tar.gz
mkdir min && cp utility/README.txt min/
for path in utility/map_*; do
  f=$(basename "$path")
  case "$f" in
    *.bz2) dec="bzcat" ; enc="bzip2 --best" ;;
    *)     dec="zcat"  ; enc="gzip -n --best" ;;
  esac
  case "$f" in
    *_uniref50.txt.gz|*_uniref90.txt.gz) $dec "$path" | head -n 1 | cut -f1-4 | $enc > "min/$f" ;;
    *)                                   $dec "$path" | head -n 1 | $enc > "min/$f" ;;
  esac
done
tar czf utility_nfDEMO.tar.gz --owner=0 --group=0 --sort=name -C min .
```

| File                       | HUMAnN option            | Contents                                                                 | Size   |
| -------------------------- | ------------------------ | ------------------------------------------------------------------------ | ------ |
| `chocophlan_nfDEMO.tar.gz` | `--nucleotide-database`  | Three ChocoPhlAn genes each of _Bacteroides dorei_ and _B. vulgatus_     | 0.9 kB |
| `uniref_nfDEMO.tar.gz`     | `--protein-database`     | UniRef90 demo proteins (FASTA and DIAMOND index)                         | 25 kB  |
| `utility_nfDEMO.tar.gz`    | utility mapping database | UniRef to KO, EC, GO, Pfam and eggNOG maps and name files, one line each | 1.8 kB |
| `demo_subset.fastq.gz`     | `--input`                | 177 HUMAnN demo reads                                                    | 6.9 kB |
| `demo_profile.txt`         | `--taxonomic-profile`    | MetaPhlAn 4 (`vJan21`) profile of the two demo species                   | 0.9 kB |

The demo databases and reads come from the HUMAnN repository (`humann/data/`, `examples/`), released
under the MIT license.
