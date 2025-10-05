# DIA-NN Test Data

This directory contains test data for the `dia_proteomics_analysis` subworkflow.

## Files

- `RD139_Narrow_UPS1_0_1fmol_inj1.mzML.tar.gz` - Compressed mzML file (66 MB)
- `REF_EColi_K12_UPS1_combined_subset_100.fasta` - E. coli proteome subset (55 KB, 100 proteins)
- `RD139_Narrow_UPS1_design.tsv` - Experimental design file

## Data Source

The original raw mass spectrometry data was obtained from:
- **URL**: https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/quantms-ci-github/MSV000087597/RD139_Narrow_UPS1_0_1fmol_inj1.raw
- **Dataset**: MSV000087597
- **Sample**: E. coli proteome spiked with UPS1 standard at 0.1 fmol
- **Method**: DIA (Data-Independent Acquisition)
- **Instrument**: Orbitrap

## Data Processing

### 1. Raw to mzML Conversion

The Thermo .raw file was converted to mzML format using ThermoRawFileParser:

```bash
ThermoRawFileParser --input RD139_Narrow_UPS1_0_1fmol_inj1.raw \
    --output . \
    --format 1
```

### 2. mzML Subsetting

To create a manageable test dataset, the full mzML file was subsetted to scans 30,000-59,999 (30,000 scans total). This range was selected because:
- It contains sufficient MS/MS data for protein identification
- It yields ~118 precursors at 1% FDR in preliminary analysis
- It produces ~10 proteins passing protein-level FDR in final quantification
- The file size is reduced from ~600 MB to ~200 MB (66 MB compressed)

Subsetting was performed using the provided Python script:

```bash
python subset_mzml.py \
    RD139_Narrow_UPS1_0_1fmol_inj1.mzML \
    RD139_Narrow_UPS1_0_1fmol_inj1_subset.mzML \
    30000 59999
```

The script uses pyOpenMS to extract the specified scan range while preserving all metadata.

### 3. Compression

The subset mzML was compressed to reduce storage and transfer time:

```bash
tar -czf RD139_Narrow_UPS1_0_1fmol_inj1.mzML.tar.gz \
    RD139_Narrow_UPS1_0_1fmol_inj1.mzML
```

**Size reduction**: 199 MB → 66 MB (67% smaller)

### 4. FASTA Database

The FASTA file contains 100 proteins from the E. coli K12 reference proteome combined with the UPS1 standard proteins. This subset was created from the full proteome to:
- Reduce library generation time during testing
- Maintain a realistic search space for FDR calculation
- Include all proteins identified in the test data

The original full database is available at:
- E. coli K12: UniProt UP000000625
- UPS1 standard: Sigma-Aldrich UPS1 (48 human proteins)

### 5. Experimental Design

The experimental design file (`RD139_Narrow_UPS1_design.tsv`) follows the quantms/DIA-NN format with two sections:

**Section 1 - Sample Information:**
```tsv
Fraction_Group	Fraction	Spectra_Filepath	Label	Sample
1	1	RD139_Narrow_UPS1_0_1fmol_inj1.mzML	1	1
```

**Section 2 - MSstats Configuration:**
```tsv
Sample	MSstats_Condition	MSstats_BioReplicate
1	UPS1_0_1fmol	1
```

Key points:
- `Spectra_Filepath` matches the mzML filename (with extension)
- Sample IDs are numeric (required by quantmsutils)
- Single-sample design for minimal test case

## Expected Test Results

With this test data, the DIA-NN analysis should produce:
- **Preliminary analysis**: ~118 precursors identified at 1% FDR
- **Individual analysis**: ~106 precursors quantified
- **Final quantification**: ~10-22 protein groups at protein-level FDR
- **MSstats output**: ~104 peptide-level quantification entries
- **mzTab output**: Protein-level quantification for identified proteins

## Test Configuration

Key parameters used for this test dataset:
- `pg_level = 1` (protein groups, not individual proteins or genes)
- `fragment_tolerance = 20 ppm`
- `precursor_tolerance = 20 ppm`
- `matrix-qvalue = 1` (permissive for test data)
- Enzyme: Trypsin
- Fixed modification: Carbamidomethyl (C)
- Variable modification: Oxidation (M)

## Notes

- This is a **minimal test dataset** designed for CI/CD and development testing
- Production analyses should use full datasets with multiple replicates
- The permissive q-value thresholds are appropriate for test validation but not for publication-quality results
- Scan range 30,000-59,999 was empirically determined to provide optimal balance between file size and identification quality

## References

- MassIVE dataset: [MSV000087597](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=MSV000087597)
- quantms pipeline: https://github.com/bigbio/quantms
- DIA-NN: https://github.com/vdemichev/DiaNN
