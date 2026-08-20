import sys
import datetime
import gzip
import pysam

def bed_to_vcf_gt(bed_file, fasta_file, vcf_file, sample_name="SAMPLE"):
    # Load the reference genome
    try:
        ref_genome = pysam.FastaFile(fasta_file)
    except OSError:
        print(f"Error: Could not open FASTA file '{fasta_file}'. Make sure it exists and is indexed (.fai).")
        sys.exit(1)

    # VCF Header
    header = [
        "##fileformat=VCFv4.2",
        f"##fileDate={datetime.date.today().strftime('%Y%m%d')}",
        f"##reference={fasta_file}",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">",
        "##ALT=<ID=BND,Description=\"Break end\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}"
    ]

    if bed_file.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt' # read text mode
    else:
        open_func = open
        mode = 'r'

    with open_func(bed_file, mode) as infile, open(vcf_file, 'w') as outfile:
        outfile.write("\n".join(header) + "\n")
        
        id_counter = 1
        
        for line in infile:
            if not line.strip() or line.startswith('#') or line.startswith('track'):
                continue
                
            parts = line.strip().split()
            if len(parts) < 3:
                continue

            chrom = parts[0]
            
            try:
                # BED is 0-based
                bed_start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            
            # VCF POS is 1-based
            pos = bed_start + 1
            
            # Fetch the ACTUAL reference base from hg38
            try:
                ref = ref_genome.fetch(chrom, bed_start, bed_start + 1).upper()
            except (KeyError, ValueError):
                # If contig not found or out of bounds, use N to prevent crash
                ref = "N"

            # Calculate SVLEN
            svlen = end - bed_start
            
            # Construct VCF fields
            vcf_id = f"bed_sv_{id_counter}"
            
            # BND Notation (REF base joined to the end position)
            alt = f"{ref}[{chrom}:{end}["
            
            # INFO field
            info = f"SVTYPE=BND;SVLEN={svlen};END={end}"
            
            # GT Logic
            format_key = "GT"
            genotype = "1/1" 
            
            outfile.write(f"{chrom}\t{pos}\t{vcf_id}\t{ref}\t{alt}\t.\tPASS\t{info}\t{format_key}\t{genotype}\n")
            
            id_counter += 1

    ref_genome.close()

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python bed2vcf.py <input.bed(.gz)> <hg38.fa> <output.vcf> [SampleName]")
    else:
        s_name = sys.argv[4] if len(sys.argv) > 4 else "SAMPLE"
        bed_to_vcf_gt(sys.argv[1], sys.argv[2], sys.argv[3], s_name)
        print(f"Conversion complete. Saved to {sys.argv[3]}")