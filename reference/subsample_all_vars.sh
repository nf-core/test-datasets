#!/usr/bin/env bash
set -euo pipefail

# Source file (change if needed)
SRC="vep_cache_113_GRCh38_chr22/homo_sapiens/113_GRCh38/22/all_vars.gz"

# Fraction to sample (1%)
FRACTION=0.001

# Output file in same directory
OUT_DIR=$(dirname "$SRC")
OUT="$OUT_DIR/subsampled_vars"

# Choose decompressor for gz or plain text
if [[ "$SRC" == *.gz ]]; then
  DECOMP="gunzip -c"
else
  DECOMP="cat"
fi

# Count lines
TOTAL_LINES=$($DECOMP "$SRC" | wc -l)
if [ "$TOTAL_LINES" -eq 0 ]; then
  echo "Input file is empty: $SRC"
  exit 1
fi

# Number of lines to sample (round to nearest integer, at least 1)
NUM_TO_SAMPLE=$(awk -v n="$TOTAL_LINES" -v f="$FRACTION" 'BEGIN{ k=int(n*f + 0.5); if(k<1) k=1; print k }')

# Temporary index file (line numbers) and cleanup
IDX_FILE=$(mktemp)
trap 'rm -f "$IDX_FILE"' EXIT

# Pick NUM_TO_SAMPLE distinct random line numbers from 1..TOTAL_LINES, then sort to preserve original order
shuf -i 1-"$TOTAL_LINES" -n "$NUM_TO_SAMPLE" | sort -n > "$IDX_FILE"

# Use awk to print only those line numbers from the original (keeps original ordering)
# Note: awk reads the index file first (NR==FNR block), then reads stdin via '-' (piped decompressed content)
$DECOMP "$SRC" | awk 'NR==FNR{keep[$1]=1; next} FNR in keep{print}' "$IDX_FILE" - > "$OUT"

echo "Wrote $NUM_TO_SAMPLE lines (approx. $(awk -v n="$TOTAL_LINES" -v k="$NUM_TO_SAMPLE" 'BEGIN{printf \"%.2f\", (k/n)*100}'))% of $TOTAL_LINES lines to: $OUT"