#!/usr/bin/env python3
"""
Subset an mzML file by scan range using pyopenms.
Usage: python subset_mzml.py input.mzML output.mzML [start_scan] [end_scan]
If only one number is provided, it takes the first N scans.
If two numbers are provided, it takes scans from start_scan to end_scan (inclusive).
"""
import sys
from pyopenms import MSExperiment, MzMLFile

def subset_mzml(input_file, output_file, start_scan=0, end_scan=None):
    """
    Subset mzML file to specified scan range.
    
    Args:
        input_file: Input mzML file path
        output_file: Output mzML file path
        start_scan: Starting scan index (0-based, inclusive)
        end_scan: Ending scan index (0-based, inclusive). If None, takes from start_scan to end of file.
    """
    print(f"Reading {input_file}...")
    
    # Load the experiment
    exp = MSExperiment()
    MzMLFile().load(input_file, exp)
    
    total_spectra = exp.getNrSpectra()
    print(f"Total spectra in file: {total_spectra}")
    
    # Determine scan range
    if end_scan is None:
        end_scan = total_spectra - 1
    
    # Validate range
    start_scan = max(0, start_scan)
    end_scan = min(end_scan, total_spectra - 1)
    
    if start_scan > end_scan:
        print(f"Error: start_scan ({start_scan}) > end_scan ({end_scan})")
        sys.exit(1)
    
    print(f"Extracting scans {start_scan} to {end_scan} ({end_scan - start_scan + 1} total)...")
    
    # Create subset
    subset_exp = MSExperiment()
    count = 0
    for i in range(start_scan, end_scan + 1):
        subset_exp.addSpectrum(exp.getSpectrum(i))
        count += 1
        if count % 1000 == 0:
            print(f"  Processed {count} scans...")
    
    # Save subset
    MzMLFile().store(output_file, subset_exp)
    print(f"✓ Wrote {count} scans (from scan {start_scan} to {end_scan}) to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python subset_mzml.py input.mzML output.mzML [start_scan] [end_scan]")
        print("  If only one number: takes first N scans")
        print("  If two numbers: takes scans from start_scan to end_scan (0-based, inclusive)")
        print("Examples:")
        print("  python subset_mzml.py input.mzML output.mzML 1000       # First 1000 scans")
        print("  python subset_mzml.py input.mzML output.mzML 1000 2000  # Scans 1000-2000")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    if len(sys.argv) == 4:
        # Single argument: treat as max_scans (backwards compatible)
        max_scans = int(sys.argv[3])
        subset_mzml(input_file, output_file, start_scan=0, end_scan=max_scans - 1)
    elif len(sys.argv) >= 5:
        # Two arguments: treat as start and end
        start_scan = int(sys.argv[3])
        end_scan = int(sys.argv[4])
        subset_mzml(input_file, output_file, start_scan=start_scan, end_scan=end_scan)
    else:
        # No arguments: default to first 20000
        subset_mzml(input_file, output_file, start_scan=0, end_scan=19999)
