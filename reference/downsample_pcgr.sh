#!/bin/bash

# PCGR Database Downsampling Script
# Filters data to chromosome 22 where possible, then samples 0.05% of entries

set -euo pipefail

SRC="20250314/data/grch38"
DST="chr22/20250314/data/grch38"

# Percentage to keep for downsampling (0.05% = 1 in 2000)
SAMPLE_RATE=0.05

echo "=== PCGR Database Downsampling Script ==="
echo "Source: $SRC"
echo "Target: $DST"
echo ""

# Check for required tools
echo "Checking required tools..."
for tool in bcftools bedtools tabix bgzip gzip awk grep zcat samtools; do
    if ! command -v "$tool" &> /dev/null; then
        echo "ERROR: Required tool '$tool' not found. Please install it first."
        exit 1
    fi
done
echo "All required tools found."
echo ""

# Create chromosome 22 region files for bedtools (separate files for different naming conventions)
CHR22_REGION_NOCHR=$(mktemp --suffix=.bed)
CHR22_REGION_CHR=$(mktemp --suffix=.bed)
echo -e "22\t0\t999999999" > "$CHR22_REGION_NOCHR"
echo -e "chr22\t0\t999999999" > "$CHR22_REGION_CHR"

# Create base directory structure
echo "Creating directory structure..."
mkdir -p "$DST"

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Filter VCF to chromosome 22 only using bcftools, then sample N%
filter_vcf_chr22() {
    local src_file="$1"
    local dst_file="$2"
    local dst_dir=$(dirname "$dst_file")
    
    mkdir -p "$dst_dir"
    echo "  Filtering VCF to chr22 + sampling ${SAMPLE_RATE}%: $(basename "$src_file")"
    
    # Try chr22 first, then 22 (different VCF naming conventions)
    # bcftools view handles both compressed and uncompressed VCFs
    # Then sample SAMPLE_RATE% of variants
    if bcftools view -r chr22 "$src_file" 2>/dev/null | grep -v "^#" | head -1 | grep -q .; then
        bcftools view -r chr22 "$src_file" | \
            awk -v pct="$SAMPLE_RATE" 'BEGIN{srand(42)} /^#/{print; next} rand()*100 < pct {print}' | \
            bgzip -c > "$dst_file"
    else
        bcftools view -r 22 "$src_file" | \
            awk -v pct="$SAMPLE_RATE" 'BEGIN{srand(42)} /^#/{print; next} rand()*100 < pct {print}' | \
            bgzip -c > "$dst_file"
    fi
    
    # Create tabix index
    tabix -p vcf "$dst_file" 2>/dev/null || true
}

# Filter BED to chromosome 22 only using bedtools, then sample N%
filter_bed_chr22() {
    local src_file="$1"
    local dst_file="$2"
    local dst_dir=$(dirname "$dst_file")
    
    mkdir -p "$dst_dir"
    echo "  Filtering BED to chr22 + sampling ${SAMPLE_RATE}%: $(basename "$src_file")"
    
    # Check which chromosome naming convention the file uses
    local first_chrom=$(zcat "$src_file" 2>/dev/null | head -1 | cut -f1)
    
    if [[ "$first_chrom" == chr* ]]; then
        # Uses "chr22" convention
        bedtools intersect -a "$src_file" -b "$CHR22_REGION_CHR" | \
            awk -v pct="$SAMPLE_RATE" 'BEGIN{srand(42)} rand()*100 < pct {print}' | \
            bgzip -c > "$dst_file"
    else
        # Uses "22" convention
        bedtools intersect -a "$src_file" -b "$CHR22_REGION_NOCHR" | \
            awk -v pct="$SAMPLE_RATE" 'BEGIN{srand(42)} rand()*100 < pct {print}' | \
            bgzip -c > "$dst_file"
    fi
    
    # Create tabix index
    tabix -p bed "$dst_file" 2>/dev/null || true
}

# Filter TSV with chromosome column to chr22 ONLY (no sampling) - for stats files
filter_tsv_chr22_nosample() {
    local src_file="$1"
    local dst_file="$2"
    local chrom_col="$3"  # Column name containing chromosome
    local dst_dir=$(dirname "$dst_file")
    
    mkdir -p "$dst_dir"
    echo "  Filtering TSV to chr22 (column: $chrom_col): $(basename "$src_file")"
    
    if [[ "$src_file" == *.gz ]]; then
        # Find column index and filter
        zcat "$src_file" | awk -F'\t' -v col="$chrom_col" '
            NR==1 {
                for(i=1; i<=NF; i++) {
                    if($i == col) { idx=i; break }
                }
                print; next
            }
            $idx == "22" || $idx == "chr22" { print }
        ' | gzip -c > "$dst_file"
    else
        awk -F'\t' -v col="$chrom_col" '
            NR==1 {
                for(i=1; i<=NF; i++) {
                    if($i == col) { idx=i; break }
                }
                print; next
            }
            $idx == "22" || $idx == "chr22" { print }
        ' "$src_file" > "$dst_file"
    fi
}

# Filter TSV with chromosome column to chr22, then sample N%
filter_tsv_chr22() {
    local src_file="$1"
    local dst_file="$2"
    local chrom_col="$3"  # Column name containing chromosome
    local dst_dir=$(dirname "$dst_file")
    
    mkdir -p "$dst_dir"
    echo "  Filtering TSV to chr22 + sampling ${SAMPLE_RATE}% (column: $chrom_col): $(basename "$src_file")"
    
    if [[ "$src_file" == *.gz ]]; then
        # Find column index, filter to chr22, then sample
        zcat "$src_file" | awk -F'\t' -v col="$chrom_col" -v pct="$SAMPLE_RATE" '
            BEGIN { srand(42) }
            NR==1 {
                for(i=1; i<=NF; i++) {
                    if($i == col) { idx=i; break }
                }
                print; next
            }
            ($idx == "22" || $idx == "chr22") && rand()*100 < pct { print }
        ' | gzip -c > "$dst_file"
    else
        awk -F'\t' -v col="$chrom_col" -v pct="$SAMPLE_RATE" '
            BEGIN { srand(42) }
            NR==1 {
                for(i=1; i<=NF; i++) {
                    if($i == col) { idx=i; break }
                }
                print; next
            }
            ($idx == "22" || $idx == "chr22") && rand()*100 < pct { print }
        ' "$src_file" > "$dst_file"
    fi
}

# Sample N% of rows from TSV (keeping header)
sample_tsv() {
    local src_file="$1"
    local dst_file="$2"
    local pct="$3"
    local dst_dir=$(dirname "$dst_file")
    
    mkdir -p "$dst_dir"
    echo "  Sampling ${pct}% of: $(basename "$src_file")"
    
    if [[ "$src_file" == *.gz ]]; then
        zcat "$src_file" | awk -v pct="$pct" '
            BEGIN { srand(42) }
            NR==1 { print; next }
            rand()*100 < pct { print }
        ' | gzip -c > "$dst_file"
    else
        awk -v pct="$pct" '
            BEGIN { srand(42) }
            NR==1 { print; next }
            rand()*100 < pct { print }
        ' "$src_file" > "$dst_file"
    fi
}

# Copy file as-is (for small reference files)
copy_file() {
    local src_file="$1"
    local dst_file="$2"
    local dst_dir=$(dirname "$dst_file")
    
    mkdir -p "$dst_dir"
    echo "  Copying: $(basename "$src_file")"
    cp "$src_file" "$dst_file"
}

# Copy tabix index if it exists
copy_index() {
    local src_file="$1"
    local dst_file="$2"
    
    if [[ -f "${src_file}.tbi" ]]; then
        cp "${src_file}.tbi" "${dst_file}.tbi"
    fi
}

# Get chr22 genes from gene_transcript_xref
get_chr22_genes() {
    zcat "$SRC/gene/tsv/gene_transcript_xref/gene_transcript_xref.tsv.gz" | \
        awk -F'\t' 'NR==1 { for(i=1;i<=NF;i++) if($i=="chrom") ci=i; if($i=="entrezgene") ei=i; if($i=="symbol") si=i; next }
                    $ci == "22" || $ci == "chr22" { print $ei"\t"$si }' | \
        sort -u
}

# Filter TSV by gene symbol (for gene-based files)
filter_by_gene_symbol() {
    local src_file="$1"
    local dst_file="$2"
    local gene_file="$3"
    local symbol_col="$4"
    local dst_dir=$(dirname "$dst_file")
    
    mkdir -p "$dst_dir"
    echo "  Filtering by chr22 genes (column: $symbol_col): $(basename "$src_file")"
    
    if [[ "$src_file" == *.gz ]]; then
        zcat "$src_file" | awk -F'\t' -v col="$symbol_col" -v gf="$gene_file" '
            BEGIN {
                while((getline line < gf) > 0) {
                    split(line, a, "\t")
                    genes[a[2]] = 1
                }
            }
            NR==1 {
                for(i=1; i<=NF; i++) {
                    if($i == col) { idx=i; break }
                }
                print; next
            }
            $idx in genes { print }
        ' | gzip -c > "$dst_file"
    else
        awk -F'\t' -v col="$symbol_col" -v gf="$gene_file" '
            BEGIN {
                while((getline line < gf) > 0) {
                    split(line, a, "\t")
                    genes[a[2]] = 1
                }
            }
            NR==1 {
                for(i=1; i<=NF; i++) {
                    if($i == col) { idx=i; break }
                }
                print; next
            }
            $idx in genes { print }
        ' "$src_file" > "$dst_file"
    fi
}

# Filter TSV by entrezgene ID
filter_by_entrezgene() {
    local src_file="$1"
    local dst_file="$2"
    local gene_file="$3"
    local entrez_col="$4"
    local dst_dir=$(dirname "$dst_file")
    
    mkdir -p "$dst_dir"
    echo "  Filtering by chr22 entrezgene (column: $entrez_col): $(basename "$src_file")"
    
    if [[ "$src_file" == *.gz ]]; then
        zcat "$src_file" | awk -F'\t' -v col="$entrez_col" -v gf="$gene_file" '
            BEGIN {
                while((getline line < gf) > 0) {
                    split(line, a, "\t")
                    if(a[1] != "") genes[a[1]] = 1
                }
            }
            NR==1 {
                for(i=1; i<=NF; i++) {
                    if($i == col) { idx=i; break }
                }
                print; next
            }
            $idx in genes { print }
        ' | gzip -c > "$dst_file"
    else
        awk -F'\t' -v col="$entrez_col" -v gf="$gene_file" '
            BEGIN {
                while((getline line < gf) > 0) {
                    split(line, a, "\t")
                    if(a[1] != "") genes[a[1]] = 1
                }
            }
            NR==1 {
                for(i=1; i<=NF; i++) {
                    if($i == col) { idx=i; break }
                }
                print; next
            }
            $idx in genes { print }
        ' "$src_file" > "$dst_file"
    fi
}

# ============================================================================
# MAIN PROCESSING
# ============================================================================

# First, extract chr22 genes for filtering gene-based files
echo ""
echo "=== Extracting chr22 gene list ==="
CHR22_GENES=$(mktemp)
get_chr22_genes > "$CHR22_GENES"
echo "Found $(wc -l < "$CHR22_GENES") genes on chr22"

# ----------------------------------------------------------------------------
# 1. ROOT FILES (copy as-is - small reference files)
# ----------------------------------------------------------------------------
echo ""
echo "=== Processing root files ==="
copy_file "$SRC/.PCGR_BUNDLE_VERSION" "$DST/.PCGR_BUNDLE_VERSION"
copy_file "$SRC/chromsize.grch38.tsv" "$DST/chromsize.grch38.tsv"
# Skip data_overview.grch38.html - it's 23MB documentation, not needed for CI
echo "  Skipping data_overview.grch38.html (large documentation file)"
copy_file "$SRC/vcf_infotags_other.tsv" "$DST/vcf_infotags_other.tsv"
copy_file "$SRC/vcf_infotags_vep.tsv" "$DST/vcf_infotags_vep.tsv"

# ----------------------------------------------------------------------------
# 2. METADATA FILES (copy as-is - small reference files)
# ----------------------------------------------------------------------------
echo ""
echo "=== Processing .METADATA files ==="
mkdir -p "$DST/.METADATA/tsv"
for f in "$SRC"/.METADATA/tsv/*.tsv; do
    copy_file "$f" "$DST/.METADATA/tsv/$(basename "$f")"
done

# ----------------------------------------------------------------------------
# 3. VCF FILES (filter to chr22)
# ----------------------------------------------------------------------------
echo ""
echo "=== Processing VCF files (filtering to chr22) ==="

for vcf_dir in "$SRC"/variant/vcf/*/; do
    dir_name=$(basename "$vcf_dir")
    for vcf_file in "$vcf_dir"*.vcf.gz; do
        if [[ -f "$vcf_file" ]]; then
            filter_vcf_chr22 "$vcf_file" "$DST/variant/vcf/$dir_name/$(basename "$vcf_file")"
        fi
    done
    # Filter varstats files to chr22 row only (they have per-chromosome stats, no sampling)
    for tsv_file in "$vcf_dir"*_varstats.tsv; do
        if [[ -f "$tsv_file" ]]; then
            filter_tsv_chr22_nosample "$tsv_file" "$DST/variant/vcf/$dir_name/$(basename "$tsv_file")" "Chromosome"
        fi
    done
done

# ----------------------------------------------------------------------------
# 4. BED FILES (filter to chr22)
# ----------------------------------------------------------------------------
echo ""
echo "=== Processing BED files (filtering to chr22) ==="

# misc/bed files
for bed_dir in "$SRC"/misc/bed/*/; do
    dir_name=$(basename "$bed_dir")
    for bed_file in "$bed_dir"*.bed.gz; do
        if [[ -f "$bed_file" ]]; then
            filter_bed_chr22 "$bed_file" "$DST/misc/bed/$dir_name/$(basename "$bed_file")"
        fi
    done
done

# gene/bed/gene_transcript_xref
mkdir -p "$DST/gene/bed/gene_transcript_xref"
for bed_file in "$SRC"/gene/bed/gene_transcript_xref/*.bed.gz; do
    if [[ -f "$bed_file" ]]; then
        filter_bed_chr22 "$bed_file" "$DST/gene/bed/gene_transcript_xref/$(basename "$bed_file")"
    fi
done
# Copy vcf_info_tags file
if [[ -f "$SRC/gene/bed/gene_transcript_xref/gene_transcript_xref.vcfanno.vcf_info_tags.txt" ]]; then
    copy_file "$SRC/gene/bed/gene_transcript_xref/gene_transcript_xref.vcfanno.vcf_info_tags.txt" \
              "$DST/gene/bed/gene_transcript_xref/gene_transcript_xref.vcfanno.vcf_info_tags.txt"
fi

# gene/bed/gene_virtual_panel - filter each panel
mkdir -p "$DST/gene/bed/gene_virtual_panel"
for bed_file in "$SRC"/gene/bed/gene_virtual_panel/*.bed.gz; do
    if [[ -f "$bed_file" ]]; then
        filter_bed_chr22 "$bed_file" "$DST/gene/bed/gene_virtual_panel/$(basename "$bed_file")"
    fi
done

# ----------------------------------------------------------------------------
# 5. VARIANT TSV FILES (filter by chr22 where possible)
# ----------------------------------------------------------------------------
echo ""
echo "=== Processing variant TSV files ==="

# clinvar.tsv.gz - has chrom column
filter_tsv_chr22 "$SRC/variant/tsv/clinvar/clinvar.tsv.gz" \
                 "$DST/variant/tsv/clinvar/clinvar.tsv.gz" "chrom"

# clinvar_sites.tsv.gz - has var_id with chr embedded (e.g., 22_12345_A_G)
mkdir -p "$DST/variant/tsv/clinvar"
echo "  Filtering clinvar_sites by chr22 var_id pattern"
zcat "$SRC/variant/tsv/clinvar/clinvar_sites.tsv.gz" | awk -F'\t' '
    NR==1 { for(i=1;i<=NF;i++) if($i=="var_id") idx=i; print; next }
    $idx ~ /^22_/ { print }
' | gzip -c > "$DST/variant/tsv/clinvar/clinvar_sites.tsv.gz"

# clinvar_oncogenic.tsv.gz - has var_id with chr embedded
echo "  Filtering clinvar_oncogenic by chr22 var_id pattern"
zcat "$SRC/variant/tsv/clinvar/clinvar_oncogenic.tsv.gz" | awk -F'\t' '
    NR==1 { for(i=1;i<=NF;i++) if($i=="var_id") idx=i; print; next }
    $idx ~ /^22_/ { print }
' | gzip -c > "$DST/variant/tsv/clinvar/clinvar_oncogenic.tsv.gz"

# clinvar_gene_varstats.tsv.gz - filter by chr22 genes
filter_by_gene_symbol "$SRC/variant/tsv/clinvar/clinvar_gene_varstats.tsv.gz" \
                      "$DST/variant/tsv/clinvar/clinvar_gene_varstats.tsv.gz" \
                      "$CHR22_GENES" "symbol"

# gwas.tsv.gz - has chromosome column
filter_tsv_chr22 "$SRC/variant/tsv/gwas/gwas.tsv.gz" \
                 "$DST/variant/tsv/gwas/gwas.tsv.gz" "chromosome"

# ----------------------------------------------------------------------------
# 6. GENE TSV FILES
# ----------------------------------------------------------------------------
echo ""
echo "=== Processing gene TSV files ==="

# gene_transcript_xref.tsv.gz - has chrom column
filter_tsv_chr22 "$SRC/gene/tsv/gene_transcript_xref/gene_transcript_xref.tsv.gz" \
                 "$DST/gene/tsv/gene_transcript_xref/gene_transcript_xref.tsv.gz" "chrom"

# gene_index.tsv.gz - filter by chr22 entrezgene
filter_by_entrezgene "$SRC/gene/tsv/gene_transcript_xref/gene_index.tsv.gz" \
                     "$DST/gene/tsv/gene_transcript_xref/gene_index.tsv.gz" \
                     "$CHR22_GENES" "ENTREZGENE"

# gene_transcript_xref_bedmap.tsv.gz - small, copy as-is
copy_file "$SRC/gene/tsv/gene_transcript_xref/gene_transcript_xref_bedmap.tsv.gz" \
          "$DST/gene/tsv/gene_transcript_xref/gene_transcript_xref_bedmap.tsv.gz"

# otp_rank.tsv.gz - filter by chr22 genes
filter_by_gene_symbol "$SRC/gene/tsv/gene_transcript_xref/otp_rank.tsv.gz" \
                      "$DST/gene/tsv/gene_transcript_xref/otp_rank.tsv.gz" \
                      "$CHR22_GENES" "symbol"

# gene_cpg - filter by chr22 genes
filter_by_gene_symbol "$SRC/gene/tsv/gene_cpg/gene_cpg.tsv.gz" \
                      "$DST/gene/tsv/gene_cpg/gene_cpg.tsv.gz" \
                      "$CHR22_GENES" "symbol"

# gene_virtual_panel - filter by chr22 genes
filter_by_gene_symbol "$SRC/gene/tsv/gene_virtual_panel/gene_virtual_panel.tsv.gz" \
                      "$DST/gene/tsv/gene_virtual_panel/gene_virtual_panel.tsv.gz" \
                      "$CHR22_GENES" "symbol"

# ----------------------------------------------------------------------------
# 7. MISC TSV FILES
# ----------------------------------------------------------------------------
echo ""
echo "=== Processing misc TSV files ==="

# cytoband.tsv.gz - has chrom column
filter_tsv_chr22 "$SRC/misc/tsv/cytoband/cytoband.tsv.gz" \
                 "$DST/misc/tsv/cytoband/cytoband.tsv.gz" "chrom"

# grantham.tsv - amino acid reference, copy as-is (small)
copy_file "$SRC/misc/tsv/grantham/grantham.tsv" "$DST/misc/tsv/grantham/grantham.tsv"

# oncogenicity.tsv - scoring rules, copy as-is (small)
copy_file "$SRC/misc/tsv/oncogenicity/oncogenicity.tsv" "$DST/misc/tsv/oncogenicity/oncogenicity.tsv"

# mutsplicedb.tsv - filter by chr22 genes
filter_by_gene_symbol "$SRC/misc/tsv/mutsplicedb/mutsplicedb.tsv" \
                      "$DST/misc/tsv/mutsplicedb/mutsplicedb.tsv" \
                      "$CHR22_GENES" "SYMBOL"

# hotspot files - filter by chr22 genes
filter_by_gene_symbol "$SRC/misc/tsv/hotspot/hotspot.tsv.gz" \
                      "$DST/misc/tsv/hotspot/hotspot.tsv.gz" \
                      "$CHR22_GENES" "symbol"
filter_by_gene_symbol "$SRC/misc/tsv/hotspot/hotspot_long.tsv.gz" \
                      "$DST/misc/tsv/hotspot/hotspot_long.tsv.gz" \
                      "$CHR22_GENES" "symbol"

# mutational_signature.tsv.gz - reference data, copy as-is
copy_file "$SRC/misc/tsv/mutational_signature/mutational_signature.tsv.gz" \
          "$DST/misc/tsv/mutational_signature/mutational_signature.tsv.gz"

# pathway.tsv.gz - filter by chr22 genes
filter_by_gene_symbol "$SRC/misc/tsv/pathway/pathway.tsv.gz" \
                      "$DST/misc/tsv/pathway/pathway.tsv.gz" \
                      "$CHR22_GENES" "symbol"

# protein_domain.tsv.gz - reference data, copy as-is (no gene column)
copy_file "$SRC/misc/tsv/protein_domain/protein_domain.tsv.gz" \
          "$DST/misc/tsv/protein_domain/protein_domain.tsv.gz"

# tmb.tsv.gz - sample data, sample 10%
sample_tsv "$SRC/misc/tsv/tmb/tmb.tsv.gz" "$DST/misc/tsv/tmb/tmb.tsv.gz" "$SAMPLE_RATE"

# Create empty clinical_trial directory (source is empty)
mkdir -p "$DST/misc/tsv/clinical_trial"

# ----------------------------------------------------------------------------
# 8. BIOMARKER TSV FILES
# ----------------------------------------------------------------------------
echo ""
echo "=== Processing biomarker TSV files ==="

mkdir -p "$DST/biomarker/tsv"

# Filter by chr22 genes (symbol column)
for f in civic.variant.tsv.gz cgi.variant.tsv.gz; do
    if [[ -f "$SRC/biomarker/tsv/$f" ]]; then
        filter_by_gene_symbol "$SRC/biomarker/tsv/$f" "$DST/biomarker/tsv/$f" "$CHR22_GENES" "symbol"
    fi
done

# Clinical and literature files - sample 10% (complex relationships)
for f in civic.clinical.tsv.gz civic.literature.tsv.gz cgi.clinical.tsv.gz cgi.literature.tsv.gz; do
    if [[ -f "$SRC/biomarker/tsv/$f" ]]; then
        sample_tsv "$SRC/biomarker/tsv/$f" "$DST/biomarker/tsv/$f" "$SAMPLE_RATE"
    fi
done

# ----------------------------------------------------------------------------
# 9. DRUG TSV FILES
# ----------------------------------------------------------------------------
echo ""
echo "=== Processing drug TSV files ==="

mkdir -p "$DST/drug/tsv"

# Drug files - sample 10%
sample_tsv "$SRC/drug/tsv/drug_all.tsv.gz" "$DST/drug/tsv/drug_all.tsv.gz" "$SAMPLE_RATE"
sample_tsv "$SRC/drug/tsv/drug_targeted.tsv.gz" "$DST/drug/tsv/drug_targeted.tsv.gz" "$SAMPLE_RATE"

# ----------------------------------------------------------------------------
# 10. PHENOTYPE TSV FILES
# ----------------------------------------------------------------------------
echo ""
echo "=== Processing phenotype TSV files ==="

mkdir -p "$DST/phenotype/tsv"

# Sample 10% of phenotype files
sample_tsv "$SRC/phenotype/tsv/phenotype_onco.tsv.gz" "$DST/phenotype/tsv/phenotype_onco.tsv.gz" "$SAMPLE_RATE"
sample_tsv "$SRC/phenotype/tsv/phenotype_umls.tsv.gz" "$DST/phenotype/tsv/phenotype_umls.tsv.gz" "$SAMPLE_RATE"

# ----------------------------------------------------------------------------
# 11. FUSION TSV FILES
# ----------------------------------------------------------------------------
echo ""
echo "=== Processing fusion TSV files ==="

mkdir -p "$DST/fusion/tsv"

# Sample 10% of fusion files
for f in "$SRC"/fusion/tsv/*.tsv.gz; do
    if [[ -f "$f" ]]; then
        sample_tsv "$f" "$DST/fusion/tsv/$(basename "$f")" "$SAMPLE_RATE"
    fi
done

# ----------------------------------------------------------------------------
# 12. EXPRESSION DATA
# ----------------------------------------------------------------------------
echo ""
echo "=== Processing expression TSV files ==="

# depmap - sample 10% of samples
mkdir -p "$DST/expression/tsv/depmap"
sample_tsv "$SRC/expression/tsv/depmap/depmap_sample_metadata.tsv.gz" \
           "$DST/expression/tsv/depmap/depmap_sample_metadata.tsv.gz" "$SAMPLE_RATE"

# For TPM files, we need to keep matched samples - extract sample IDs first
echo "  Processing depmap TPM (sampling columns)..."
# Get sampled sample IDs
DEPMAP_SAMPLES=$(mktemp)
zcat "$DST/expression/tsv/depmap/depmap_sample_metadata.tsv.gz" | \
    awk -F'\t' 'NR>1 {print $1}' > "$DEPMAP_SAMPLES"

# Filter TPM columns to match sampled samples
zcat "$SRC/expression/tsv/depmap/depmap_tpm.tsv.gz" | awk -F'\t' -v sf="$DEPMAP_SAMPLES" '
    BEGIN {
        while((getline s < sf) > 0) samples[s] = 1
    }
    NR==1 {
        printf $1
        for(i=2; i<=NF; i++) {
            if($i in samples) {
                keep[i] = 1
                printf "\t%s", $i
            }
        }
        printf "\n"
        next
    }
    {
        printf $1
        for(i=2; i<=NF; i++) {
            if(i in keep) printf "\t%s", $i
        }
        printf "\n"
    }
' | gzip -c > "$DST/expression/tsv/depmap/depmap_tpm.tsv.gz"
rm "$DEPMAP_SAMPLES"

# treehouse - sample 10%
mkdir -p "$DST/expression/tsv/treehouse"
sample_tsv "$SRC/expression/tsv/treehouse/treehouse_sample_metadata.tsv.gz" \
           "$DST/expression/tsv/treehouse/treehouse_sample_metadata.tsv.gz" "$SAMPLE_RATE"

# For TPM, filter columns to match sampled samples
echo "  Processing treehouse TPM (sampling columns)..."
TREEHOUSE_SAMPLES=$(mktemp)
zcat "$DST/expression/tsv/treehouse/treehouse_sample_metadata.tsv.gz" | \
    awk -F'\t' 'NR>1 {print $1}' > "$TREEHOUSE_SAMPLES"

zcat "$SRC/expression/tsv/treehouse/treehouse_tpm.tsv.gz" | awk -F'\t' -v sf="$TREEHOUSE_SAMPLES" '
    BEGIN {
        while((getline s < sf) > 0) samples[s] = 1
    }
    NR==1 {
        printf $1
        for(i=2; i<=NF; i++) {
            if($i in samples) {
                keep[i] = 1
                printf "\t%s", $i
            }
        }
        printf "\n"
        next
    }
    {
        printf $1
        for(i=2; i<=NF; i++) {
            if(i in keep) printf "\t%s", $i
        }
        printf "\n"
    }
' | gzip -c > "$DST/expression/tsv/treehouse/treehouse_tpm.tsv.gz"
rm "$TREEHOUSE_SAMPLES"

# TCGA - process each cancer type
mkdir -p "$DST/expression/tsv/tcga"
echo "  Processing TCGA expression files..."

for metadata_file in "$SRC"/expression/tsv/tcga/*_sample_metadata.tsv*; do
    if [[ -f "$metadata_file" ]]; then
        base=$(basename "$metadata_file")
        cancer_type=$(echo "$base" | sed 's/_sample_metadata.tsv.*//' | sed 's/_sample_metadata.tsv//')
        
        # Check if file is gzipped or not
        if [[ "$metadata_file" == *.gz ]]; then
            ext=".tsv.gz"
            sample_tsv "$metadata_file" "$DST/expression/tsv/tcga/${cancer_type}_sample_metadata.tsv.gz" "$SAMPLE_RATE"
            
            # Get sample IDs
            TCGA_SAMPLES=$(mktemp)
            zcat "$DST/expression/tsv/tcga/${cancer_type}_sample_metadata.tsv.gz" | \
                awk -F'\t' 'NR>1 {print $1}' > "$TCGA_SAMPLES"
            
            # Filter TPM
            tpm_file="$SRC/expression/tsv/tcga/${cancer_type}_tpm.tsv.gz"
            if [[ -f "$tpm_file" ]]; then
                echo "    Filtering ${cancer_type}_tpm.tsv.gz"
                zcat "$tpm_file" | awk -F'\t' -v sf="$TCGA_SAMPLES" '
                    BEGIN {
                        while((getline s < sf) > 0) samples[s] = 1
                    }
                    NR==1 {
                        printf $1
                        for(i=2; i<=NF; i++) {
                            if($i in samples) {
                                keep[i] = 1
                                printf "\t%s", $i
                            }
                        }
                        printf "\n"
                        next
                    }
                    {
                        printf $1
                        for(i=2; i<=NF; i++) {
                            if(i in keep) printf "\t%s", $i
                        }
                        printf "\n"
                    }
                ' | gzip -c > "$DST/expression/tsv/tcga/${cancer_type}_tpm.tsv.gz"
            fi
            rm "$TCGA_SAMPLES"
        else
            # Uncompressed file
            sample_tsv "$metadata_file" "$DST/expression/tsv/tcga/${cancer_type}_sample_metadata.tsv" "$SAMPLE_RATE"
            
            # Get sample IDs
            TCGA_SAMPLES=$(mktemp)
            awk -F'\t' 'NR>1 {print $1}' "$DST/expression/tsv/tcga/${cancer_type}_sample_metadata.tsv" > "$TCGA_SAMPLES"
            
            # Filter TPM
            tpm_file="$SRC/expression/tsv/tcga/${cancer_type}_tpm.tsv"
            if [[ -f "$tpm_file" ]]; then
                echo "    Filtering ${cancer_type}_tpm.tsv"
                awk -F'\t' -v sf="$TCGA_SAMPLES" '
                    BEGIN {
                        while((getline s < sf) > 0) samples[s] = 1
                    }
                    NR==1 {
                        printf $1
                        for(i=2; i<=NF; i++) {
                            if($i in samples) {
                                keep[i] = 1
                                printf "\t%s", $i
                            }
                        }
                        printf "\n"
                        next
                    }
                    {
                        printf $1
                        for(i=2; i<=NF; i++) {
                            if(i in keep) printf "\t%s", $i
                        }
                        printf "\n"
                    }
                ' "$tpm_file" > "$DST/expression/tsv/tcga/${cancer_type}_tpm.tsv"
            fi
            rm "$TCGA_SAMPLES"
        fi
    fi
done

# ----------------------------------------------------------------------------
# 13. MISC OTHER FILES
# ----------------------------------------------------------------------------
echo ""
echo "=== Processing misc other files ==="

# misc/fasta - copy reference files as-is (needed for functionality)
# misc/fasta - extract chr22 only from genome assembly (full genome is 855MB!)
if [[ -d "$SRC/misc/fasta" ]]; then
    mkdir -p "$DST/misc/fasta/ancestor"
    mkdir -p "$DST/misc/fasta/assembly"
    
    # Extract chr22 from the genome assembly using samtools
    GENOME_FA="$SRC/misc/fasta/assembly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    if [[ -f "$GENOME_FA" ]]; then
        echo "  Extracting chr22 from genome assembly..."
        # Use samtools faidx to extract chr22
        samtools faidx "$GENOME_FA" 22 | bgzip -c > "$DST/misc/fasta/assembly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        # Create new index for the chr22-only fasta
        samtools faidx "$DST/misc/fasta/assembly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    fi
fi

# misc/other - copy only essential files, skip large HTML documentation
if [[ -d "$SRC/misc/other" ]]; then
    mkdir -p "$DST/misc/other/msi_classification"
    # Copy only the RDS model file, skip the HTML (1.5MB documentation)
    if [[ -f "$SRC/misc/other/msi_classification/tcga_msi_classifier.rds" ]]; then
        copy_file "$SRC/misc/other/msi_classification/tcga_msi_classifier.rds" \
                  "$DST/misc/other/msi_classification/tcga_msi_classifier.rds"
    fi
    echo "  Skipping tcga_msi_classifier.html (large documentation file)"
fi

# misc/bed/ncer - filter if exists
if [[ -d "$SRC/misc/bed/ncer" ]]; then
    mkdir -p "$DST/misc/bed/ncer"
    for bed_file in "$SRC"/misc/bed/ncer/*.bed.gz; do
        if [[ -f "$bed_file" ]]; then
            filter_bed_chr22 "$bed_file" "$DST/misc/bed/ncer/$(basename "$bed_file")"
        fi
    done
fi

# Cleanup
rm "$CHR22_GENES"
rm "$CHR22_REGION_NOCHR"
rm "$CHR22_REGION_CHR"

# ----------------------------------------------------------------------------
# SUMMARY
# ----------------------------------------------------------------------------
echo ""
echo "=== Downsampling Complete ==="
echo ""
echo "Source size:"
du -sh "$SRC"
echo ""
echo "Target size:"
du -sh "$DST"
echo ""
echo "File counts:"
echo "  Source: $(find "$SRC" -type f | wc -l) files"
echo "  Target: $(find "$DST" -type f | wc -l) files"
