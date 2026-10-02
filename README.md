# test-datasets: `plasmodiumdrugres`

This branch contains test data to be used for automated testing with the [nf-core/plasmodiumdrugres](https://github.com/nf-core/plasmodiumdrugres) pipeline.

## Content of this repository

The pipeline analyzes _Plasmodium_ drug resistance markers from allele tables or [Portable Microhaplotype Object (PMO)](https://plasmogenepi.github.io/PMO_Docs/) files, translating loci of interest to amino acid changes and estimating single-locus and multi-locus allele frequencies and prevalences.

All files live under `testdata/`. Files used by the minimal `test`, `test_pmo` and nf-test profiles are realistic but simulated data. Files prefixed with `dataset1_` are real, published data used by the `test_full` profile (see [Full-size test data](#full-size-test-data-test_full)).

### Pipeline inputs

- `testdata/example_PMO.json`: Example Portable Microhaplotype Object (PMO) file for the PMO entry point including simulated data.
- `testdata/allele_table.tsv`: Microhaplotype allele table (specimen, target, allele, reads, sequence) for the allele-table entry point extracted from `example_PMO.json` (see extract_allele_table.nf module).
- `testdata/panel_info.bed`: Panel information BED including insert reference sequences extracted from `example_PMO.json` (see extract_panel_info_to_bed.nf module).
- `testdata/panel_info_no_ref.bed`: Panel information BED without reference sequences extracted from `example_PMO.json` (see extract_panel_info_to_bed.nf module).
- `testdata/dummy_panel_info_fake_chroms.bed`: Panel information BED with synthetic chromosome names to test adding ref from whole genome using `insert_refseqs.fasta`.
- `testdata/insert_refseqs.fasta`: Targeted insert/reference sequences matching panel target names.
- `testdata/loci_of_interest.bed`: Drug-resistance loci of interest (amino acid positions) used for translation.
- `testdata/loci_groups.tsv`: Multi-locus groups (e.g. `pfdhfr_pfdhps`, `crt`) for multi-locus frequency estimation.
- `testdata/population_assignment.tsv`: Specimen-to-population assignment table (see extract_population_map_from_pmo.nf module).
- `testdata/population_assignment_unassigned.tsv`: Specimen-to-population assignment table missing some specimens (`Vietnam2018-23`, `Vietnam2018-24`) to test handling of unassigned specimens.

### Intermediate and module test inputs

- `testdata/amino_acid_calls.tsv`: Amino acid call table produced from translation of loci of interest.
- `testdata/loci_of_interest_mhaps.bed`: Loci-of-interest annotations linked to microhaplotype targets (for microhaplotype-based SLAF tests).
- `testdata/mhaps_slaf.tsv`: Microhaplotype single-locus allele frequencies.
- `testdata/aa_mlaf.tsv`: Amino acid multi-locus allele frequencies (MLAF).
- `testdata/mlaf_pop1.tsv`: Multi-locus allele frequencies for population `pop1`.
- `testdata/mlaf_pop2.tsv`: Multi-locus allele frequencies for population `pop2`.
- `testdata/population_assignment_indexed.tsv`: Specimen-to-population assignment including population index IDs.
- `testdata/population_assignment_indexed_unassigned.tsv`: Specimen-to-population assignment including population index IDs, missing some specimens (`Vietnam2018-23`, `Vietnam2018-24`) to test handling of unassigned specimens.
- `testdata/population_index_lookup.tsv`: Lookup table mapping population index IDs to population names.
- `testdata/empty_population_index_lookup.tsv`: Empty population index lookup file for edge-case tests.
- `testdata/allele_prev.tsv`: Expected allele prevalence estimates.
- `testdata/aa_slaf.tsv`: Expected amino acid single-locus allele frequencies (SLAF).
- `testdata/amino_acid_calls.aa_sl_from_ml.tsv`: Expected single-locus amino acid frequencies derived from multi-locus frequencies.
- `testdata/sl_from_ml.tsv`: Expected single-locus frequencies derived from multi-locus allele frequencies.
- `testdata/sl_from_mlaf_pop1.tsv`: Expected single-locus frequencies derived from MLAF for population `pop1`.
- `testdata/sl_from_mlaf_pop2.tsv`: Expected single-locus frequencies derived from MLAF for population `pop2`.
- `testdata/sl_pop1.tsv`: Expected single-locus allele frequencies and prevalences for population `pop1`.
- `testdata/sl_pop2.tsv`: Expected single-locus allele frequencies and prevalences for population `pop2`.

### Full-size test data (`test_full`)

Real _Plasmodium falciparum_ genomic surveillance data from Eswatini, Namibia, South Africa and Zambia, generated with the MAD4HatTeR amplicon sequencing panel (Aranda-Díaz _et al._ 2025b) on specimens collected between 2022 and 2024. This is "Dataset 1" from the PMO paper, where it was used to estimate drug resistance marker prevalence by country and province:

> Hathaway NJ, Murie K, Murphy M _et al._ The Portable Microhaplotype Object and Tools. bioRxiv 2025. DOI: [10.64898/2025.12.10.693568](https://doi.org/10.64898/2025.12.10.693568)

The original PMO (`dataset1_pmo_qc_pass.json.gz`) is archived on Zenodo: [10.5281/zenodo.20550920](https://doi.org/10.5281/zenodo.20550920) (CC-BY-4.0). It contains 1,646 specimens that passed QC (Zambia 916, South Africa 476, Eswatini 144, Namibia 110), 81 targets across panels A52 and AB2, and microhaplotypes called with MAD4HatTeR pipeline v0.2.1. Only the drug resistance targets of the panel were made publicly available. The specimens come from these studies:

- Aranda-Díaz A, Mwanza S, Makhanthisa TI _et al._ _Plasmodium falciparum_ Genomic Surveillance Reveals a Diversity of Kelch 13 Mutations in Zambia. Am J Trop Med Hyg 2025a. DOI: [10.4269/ajtmh.25-0110](https://doi.org/10.4269/ajtmh.25-0110)
- Eloff L, Aranda-Díaz A, Routledge I _et al._ High Prevalence of Molecular Markers Associated with Artemisinin, Sulphadoxine and Pyrimethamine Resistance in Northern Namibia. medRxiv 2025. DOI: [10.1101/2025.01.09.25320247](https://doi.org/10.1101/2025.01.09.25320247)
- Nhlengethwa N, Aranda-Díaz A, Vilakati S _et al._ Genomic Surveillance Reveals Clusters of _Plasmodium falciparum_ Antimalarial Resistance Markers in Eswatini, a Low-Transmission Setting. medRxiv 2025. DOI: [10.1101/2025.07.30.25332463](https://doi.org/10.1101/2025.07.30.25332463)
- Raman J, Mabona M, Nyawo Q _et al._ Very low prevalence of validated kelch13 mutations and absence of hrp2/3 double gene deletions in South African malaria-eliminating districts (2022-2024). medRxiv 2025. DOI: [10.1101/2025.03.31.25324948](https://doi.org/10.1101/2025.03.31.25324948)

Panel reference: Aranda-Díaz A, Neubauer Vickers E, Murie K _et al._ Sensitive and modular amplicon sequencing of _Plasmodium falciparum_ diversity and resistance for research and public health. Sci Rep 2025b;15:10737. DOI: [10.1038/s41598-025-94716-5](https://doi.org/10.1038/s41598-025-94716-5)

Files:

- `testdata/dataset1_pmo_qc_pass_with_locations.json.gz`: The Zenodo `dataset1_pmo_qc_pass.json.gz` with an `insert_location` (chrom, 0-based start, end, strand, `genome_id` 0 = 3D7 PlasmoDB v65, and insert `ref_seq`) added to every `target_info` entry. The pipeline needs these locations to build the panel BED. Locations and reference sequences come from the MAD4HatTeR panel insert BED; for all 81 targets the `ref_seq` matches an observed microhaplotype exactly. The order of `panel_info` is also swapped so that AB2 (all 81 targets) comes before A52 (59 targets), with `library_sample_info[].panel_id` remapped to match. This works around pmotools-python 1.1.0 `extract_insert_of_panels` exporting only the first panel. Nothing else in the PMO is changed.
- `testdata/dataset1_loci_of_interest.bed`: The loci reported in the paper, _crt_ (`PF3D7_0709000`) codons 76 and 97 and _k13_ (`PF3D7_1343700`) codon 441, plus _dhfr_ (`PF3D7_0417200`) codons 51, 59 and 108 for multi-locus estimation.
- `testdata/dataset1_loci_groups.tsv`: Groups _dhfr_ codons 51, 59 and 108 (the triple mutant haplotype) for multi-locus estimation.
