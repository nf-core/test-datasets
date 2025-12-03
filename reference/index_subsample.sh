#!/usr/bin/env bash
set -euo pipefail

SRC="${1:-vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/22/subsampled_vars}"
OUT="${2:-${SRC%.gz}.sorted.bgz}"   # default: same name without .gz + .sorted.bgz

# temp dir for sort (change if you want)
TMPDIR="${TMPDIR:-/tmp}"

if [ ! -r "$SRC" ]; then
  echo "Input file not found or not readable: $SRC" >&2
  exit 2
fi

# choose decompressor (works for .gz or plain files)
if [[ "$SRC" == *.gz ]]; then
  DECOMP="gunzip -c"
else
  DECOMP="cat"
fi

echo "Reading: $SRC"
echo "Writing compressed sorted output: $OUT"
echo "Index will be created next to it: ${OUT}.csi"

# Pipeline:
# 1) Decompress (if gz)
# 2) Ensure column 6 has a numeric end (set to start col5 if '.' or empty)
# 3) Sort by col1 (chr) then numeric col5 (start)
# 4) bgzip write to $OUT
# 5) tabix index with chr=1 start=5 end=6
$DECOMP "$SRC" \
  | awk -F'\t' 'BEGIN{OFS=FS} { if($6=="" || $6==".") $6=$5; print }' \
  | LC_ALL=C sort -k1,1 -k5,5n -T "$TMPDIR" \
  | bgzip -c > "$OUT"

# create (or overwrite) CSI index; chr column 1, start 5, end 6
# Use --csi to request a CSI index (output will be ${OUT}.csi)
tabix -f --csi -s 1 -b 5 -e 6 "$OUT"

echo "Done: $OUT and ${OUT}.csi"