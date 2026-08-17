# test-datasets: `raredisease`

This branch contains test data to be used for automated testing with the [nf-core/raredisease](https://github.com/nf-core/raredisease) pipeline.

## Content of this repository

`reference/`: background resources needed by tools of raredisease pipeline

`testdata/`: chr20 test resources

`reference/grch38_gnomad_reformated_-r3.1.1-.vcf.gz`: Gnomad vcf file containing entries for the region chr20:90000-92000

`reference/grch38_vcfanno_config_-v0.2-_chr20.toml`: TOML file for small test

`reference/vcfanno_grch38_small_test.tar.gz`: the archived files of grch38_*.{vcf.gz, vcf.gz.tbi} for small test

`reference/genome.ploidy_priors.tsv`: Contains contig ploidy priors for gatk4's DetermineGermlineContigPloidy

`reference/genome.ploidy_model.tar.gz`: tar gzipped directory containing the ploidy model files

`reference/genome.germline_cnv_model.tar.gz`: tar gzipped directory containing the cnv model files

`reference/mobile_elemement_references.tsv`: tsv file with paths to the mobile element locations on chromosome 21

### For SNV subworkflow, scoring variants with MIVMIR, GICAM
`reference/rank_model_dummy_mivmir_gicam_unit_test.ini`: Dummy genmod config for running genmod scoring

`testdata/justhusky_snv_rank_variants_mivmir_gicam.vcf`: Data for running `rank_variants` subworkflow with MIVMIR, GICAM

`testdata/justhusky_snv_mivmir.vcf`: Data for running MIVMIR module test

`testdata/justhusky_snv_gicam.vcf`: Data for running GICAM module test

### For Mitochondrial subworkflow

`reference/Homo_sapiens_assembly38_chr20_chrM.fasta`: chr20 and chrM hg38 reference fasta file

`reference/Homo_sapiens_assembly38_chr20_chrM.fasta.fai`: chr20 and chrM hg38 reference index fasta file

`reference/hg38.chrM.fa`: chrM hg38 reference fasta file

`reference/hg38.chrM.fa.fai`: chrM hg38 reference index fasta file

`reference/hg38.chrMshifted8000.fa`: chrM hg38 reference fasta file shifted by 8000 bp

`reference/hg38.chrMshifted8000.fa.fai`: chrM hg38 reference index fasta file shifted by 8000 bp

`reference/control_region_shifted.chrM.interval_list`: Subset of mitochondrial control regions shifted by 8000 bp

`reference/non_control_region.chrM.interval_list`: Subset of mitochondrial non control regions

`reference/ShiftBack.chain`: Used to liftover from shifted regions

`testdata/NA12878_mito_1.fq.gz`: Test fastq file 1 with chr20 and chrM reads

`testdata/NA12878_mito_2.fq.gz`: Test fastq file 2 with chr20 and chrM reads

`testdata/NA12878_sorted_chrM_chr20_rehead_60pdown.cram`: Test alignment file with chrM and chr20 downsampled and compressed 

`testdata/samplesheet_MT.csv`: samplesheet containing fastq from chr20 and chrM

`testdata/samplesheet_2_samples.csv`: samplesheet containing fastq from chr20 and chrM for 2 patients

## Minimal test dataset

A ~40MB, 9-locus, 3-sample trio dataset with a sliced local reference, used to drive the `test_minimal` CI profile (GTF-mode VEP, no Ensembl cache). Regions/labels are defined in `manifests/fasta_region_table.tsv` and their real-genome origin is documented in `manifests/coordinate_provenance.tsv`.

Loci covered, one per region id used throughout `fastq/`: `01_chrM` (MT_whole, MT SNVs), `02_chrX` (SOX3, STR), `03_chr12` (RILPL1, STR), `04_chr20b` (ME_ALU_L1_HERV, mobile element), `05_chr7` (ME_SVA, mobile element), `06_chr20a` (NOP56, STR), `07_chr16` (XYLT1, STR), `08_chr21` (CSTB, STR), `09_chr1` (INV_DUP, SV).

### `fastq/`

Trio of samples (`ACC13778A1`, `ACC13778A2`, `ACC13778A3`), each with the same 9 regions x R1/R2/singleton layout:

`fastq/<sample>/<region>_R1.fastq.gz`: forward reads extracted from the real-genome region for this locus/sample

`fastq/<sample>/<region>_R2.fastq.gz`: reverse reads, paired with the R1 above

`fastq/<sample>/<region>_singleton.fastq.gz`: reads whose mate fell outside the sliced region during extraction

### `alignment/`

Pre-aligned BAM/CRAM versions of the trio's `fastq/` reads, for the `test_align`/`test_align_singleton` profiles (which skip alignment and feed the pipeline pre-aligned files directly). All derived from the same underlying alignment as `fastq/`, so results should match the fastq-driven `test`/`test_singleton` profiles.

`alignment/ACC13778A1_sorted_md.cram` (+ `.crai`): proband, pre-aligned and CRAM-compressed, used by `test_align`'s trio (mixed BAM/CRAM input) to also exercise CRAM ingestion

`alignment/ACC13778A1_sorted_md.bam` (+ `.bai`): proband, pre-aligned BAM, converted from the CRAM above (same underlying alignment -- alignment is sample-level and doesn't depend on which other samples are present in a run); used as the sole input for `test_align_singleton`

`alignment/ACC13778A2_sorted_md.bam` (+ `.bai`): father, pre-aligned BAM, used by `test_align`'s trio

`alignment/ACC13778A3_sorted_md.bam` (+ `.bai`): mother, pre-aligned BAM, used by `test_align`'s trio

### `precalled/`

Case-level precalled VCFs for case `amusingmarmoset`, one per variant type, nested by scope: `trio/` for the `test_vcf` profile, `singleton/` for `test_vcf_singleton` (the two scopes' calls genuinely differ -- jointly-called across the family vs. the proband alone -- unlike `alignment/` above, which has no such fork). These profiles skip calling in favor of the precalled data. Only the proband (`ACC13778A1`) gets a samplesheet row per type -- parents (trio scope) are referenced via `paternal_id`/`maternal_id` only, per the pipeline's precalled-VCF samplesheet convention (see `docs/usage.md` in nf-core/raredisease).

`precalled/trio/amusingmarmoset_precalled_snv.vcf.gz` (+ `.tbi`): case-level SNV calls

`precalled/trio/amusingmarmoset_precalled_sv.vcf.gz` (+ `.tbi`): case-level SV calls

`precalled/trio/amusingmarmoset_precalled_mt.vcf.gz` (+ `.tbi`): case-level MT calls

`precalled/trio/amusingmarmoset_precalled_me.vcf.gz` (+ `.tbi`): case-level mobile-element calls

`precalled/trio/amusingmarmoset_precalled_repeat.vcf.gz` (+ `.tbi`): case-level repeat-expansion calls, pre-Stranger annotation (matches `CALL_REPEAT_EXPANSIONS`'s own `vcf` emit, the convention the pipeline expects for this precalled type)

`precalled/singleton/amusingmarmoset_precalled_snv.vcf.gz` (+ `.tbi`), `..._sv.vcf.gz`, `..._mt.vcf.gz`, `..._me.vcf.gz`, `..._repeat.vcf.gz` (+ `.tbi` each): same 5 types, called for the proband alone (no parents)

### `manifests/`

`manifests/samplesheet.csv`: pipeline-format samplesheet (trio, sex/phenotype/pedigree) pointing at the `fastq/` files

`manifests/samplesheet_singleton.csv`: samplesheet for the proband alone (no parents), pointing at the `fastq/` files, for the `test_singleton` profile

`manifests/samplesheet_align.csv`: trio samplesheet pointing at `alignment/` (2 BAM parents + 1 CRAM proband), for the `test_align` profile

`manifests/samplesheet_align_singleton.csv`: singleton samplesheet pointing at `alignment/`'s proband BAM, for the `test_align_singleton` profile

`manifests/samplesheet_vcf.csv`: trio samplesheet (proband-only rows) pointing at `precalled/trio/`, for the `test_vcf` profile

`manifests/samplesheet_vcf_singleton.csv`: singleton samplesheet pointing at `precalled/singleton/`, for the `test_vcf_singleton` profile

`manifests/coordinate_provenance.tsv`: maps each local sliced contig/region back to its real GRCh38 chrom:start-end origin, with variant types present and notes

`manifests/fasta_region_table.tsv`: the 9 regions' local coordinates, lengths, and locus labels used to build `minimal_reference.fasta`

`manifests/fastq_extraction_manifest.tsv`: per-region, per-sample read/singleton counts and source fastq paths from the extraction step

`manifests/mt_downsample_manifest.tsv`: before/after read count, coverage, and file size for downsampling each sample's chrM reads to ~150x

`manifests/param_to_file_mapping.tsv`: maps each pipeline config param (fasta, fai, sequence_dictionary, etc.) to its minimized file and original vs. minimized size

### `reference_sliced/`

`reference_sliced/minimal_reference.fasta`: the sliced multi-contig reference (9 regions, ~226KB), original contig names kept, local coordinates

`reference_sliced/minimal_reference.fasta.fai`: samtools faidx index for the sliced reference

`reference_sliced/minimal_reference.dict`: samtools/Picard sequence dictionary for the sliced reference

`reference_sliced/minimal_reference.fasta.amb`: BWA index, ambiguous-base positions

`reference_sliced/minimal_reference.fasta.ann`: BWA index, sequence/contig annotation

`reference_sliced/minimal_reference.fasta.pac`: BWA index, packed 2-bit sequence

`reference_sliced/minimal_reference.fasta.bwt.2bit.64`: BWA-MEM2 index, 2-bit Burrows-Wheeler transform

`reference_sliced/minimal_reference.fasta.0123`: BWA-MEM2 index, 2-bit packed sequence

### `resources_remapped/beds/`

`resources_remapped/beds/target_bed.remapped.bed`: target capture regions, remapped to local coordinates

`resources_remapped/beds/sambamba_regions.remapped.bed`: coverage-QC regions for sambamba, remapped

`resources_remapped/beds/manta_call_regions.remapped.bed.gz` (+ `.tbi`): Manta SV-calling regions, remapped and indexed

`resources_remapped/beds/intervals_wgs.interval_list`: whole-genome calling intervals for the sliced reference

`resources_remapped/beds/intervals_y.interval_list`: chrY-equivalent calling intervals for the sliced reference

`resources_remapped/beds/GRCh38_PAR.remapped.bed`: pseudoautosomal-region bed, remapped (empty for this dataset, no PAR overlap in the sliced regions)

### `resources_remapped/mobile_elements/`

`resources_remapped/mobile_elements/grch38_alu.remapped.bed`: ALU mobile-element reference positions, remapped

`resources_remapped/mobile_elements/grch38_herv.remapped.bed`: HERV mobile-element reference positions, remapped

`resources_remapped/mobile_elements/grch38_l1.remapped.bed`: L1 mobile-element reference positions, remapped

`resources_remapped/mobile_elements/grch38_sva.remapped.bed`: SVA mobile-element reference positions, remapped

`resources_remapped/mobile_elements/grch38_alu_swegen.remapped.vcf.gz`: SweGen ALU frequency annotation VCF, remapped

`resources_remapped/mobile_elements/grch38_herv_swegen.remapped.vcf.gz`: SweGen HERV frequency annotation VCF, remapped

`resources_remapped/mobile_elements/grch38_l1_swegen.remapped.vcf.gz`: SweGen L1 frequency annotation VCF, remapped

`resources_remapped/mobile_elements/grch38_sva_swegen.remapped.vcf.gz`: SweGen SVA frequency annotation VCF, remapped

`resources_remapped/mobile_elements/mobile_element_references_remapped.tsv`: maps each ME type (ALU/HERV/L1/SVA) to its remapped bed file, for `--mobile_element_references`

`resources_remapped/mobile_elements/mobile_element_svdb_annotations_remapped.csv`: maps each SweGen ME VCF to its SVDB frequency/count INFO key rename, for `--mobile_element_svdb_annotations`

### `resources_remapped/peddy/`

`resources_remapped/peddy/custom_sites.txt`: 592 chrom:pos:ref:alt sites within the sliced regions, used as `--peddy_sites` since peddy's bundled hg19/hg38 sites don't overlap the local coordinates

### `resources_remapped/rank_models/`

`resources_remapped/rank_models/grch38_rank_model_-v0.4-.minimal_trimmed.ini`: production SNV/MT rank model with 13 plugin-dependent sections removed (CSQ fields absent in GTF-mode VEP output would otherwise crash genmod scoring)

### `resources_remapped/svdb/`

`resources_remapped/svdb/grch38_gnomad.v4.1.sv.sites.no_cnv.remapped.vcf.gz` (+ `.tbi`): gnomAD v4.1 SV sites for the sliced regions, for SVDB frequency annotation

`resources_remapped/svdb/grch38_clinvar_sv_pathogenic_-20250728-.remapped.bedpe`: ClinVar pathogenic SV records overlapping the sliced regions

`resources_remapped/svdb/grch38_clinvar_sv_non_pathogenic_-20250728-.remapped.bedpe`: ClinVar non-pathogenic SV records overlapping the sliced regions

`resources_remapped/svdb/svdb_query_dbs_remapped.csv`: maps the gnomAD-SV VCF to its SVDB frequency/count INFO key rename, for `--svdb_query_dbs`

`resources_remapped/svdb/svdb_query_bedpedbs_remapped.csv`: maps the two ClinVar bedpe files to their SVDB frequency/count INFO key renames, for `--svdb_query_bedpedbs`

### `resources_remapped/variant_catalog/`

`resources_remapped/variant_catalog/grch38_expansionhunter_variant_catalog_remapped.json`: ExpansionHunter STR catalog, trimmed to the 4 STR loci (SOX3, RILPL1, NOP56, XYLT1) in local coordinates

### `resources_remapped/vcfanno/`

`resources_remapped/vcfanno/grch38_gnomad_reformatted_merged_r4.1.remapped.vcf.gz` (+ `.tbi`): gnomAD v4.1 genome AF, remapped, used as `--gnomad_af` source for vcfanno

`resources_remapped/vcfanno/grch38_gnomad_af.remapped.tab.gz` (+ `.tbi`): gnomAD AF reformatted as a tabix'd tab file for vcfanno's `--gnomad_af`

`resources_remapped/vcfanno/grch38_gnomad_genomes_mt_-r3.1-.vcf.gz` (+ `.tbi`): gnomAD MT-specific genome AF, remapped

`resources_remapped/vcfanno/grch28_grch37_genbank_haplogroup_-2015-08-01-.vcf.gz` (+ `.tbi`): GenBank mtDNA haplogroup reference VCF, remapped

`resources_remapped/vcfanno/grch38_vcfanno_config_remapped.toml`: vcfanno TOML config wiring the above resources to the sliced/local-coordinate contigs

`resources_remapped/vcfanno/vcfanno_resources_manifest.txt`: list of vcfanno resource file paths referenced by the TOML config, for `--vcfanno_resources`

### `resources_remapped/vep_gtf/`

`resources_remapped/vep_gtf/minimal_reference.sorted.gtf.gz` (+ `.tbi`): sorted, indexed GTF for the sliced reference, VEP's GTF-mode annotation source, replacing the Ensembl cache

`resources_remapped/vep_gtf/grch38_gnomad_pli_per_gene_-r4.1-.txt`: gnomAD v4.1 per-gene pLI scores, for the VEP pLI plugin

`resources_remapped/vep_gtf/pLI.pm`: VEP pLI plugin Perl module

`resources_remapped/vep_gtf/pLI_values.txt`: pLI plugin's score lookup table (subset covering genes in-scope)

`resources_remapped/vep_gtf/LoFtool.pm`: VEP LoFtool plugin Perl module

`resources_remapped/vep_gtf/LoFtool_scores.txt`: LoFtool plugin's gene-level loss-of-function intolerance scores

`resources_remapped/vep_gtf/SpliceAI.pm`: VEP SpliceAI plugin Perl module

`resources_remapped/vep_gtf/spliceai_snv.vcf.gz` (+ `.tbi`): precomputed SpliceAI SNV scores for the sliced regions

`resources_remapped/vep_gtf/spliceai_indel.vcf.gz` (+ `.tbi`): precomputed SpliceAI indel scores for the sliced regions

`resources_remapped/vep_gtf/plugin_config.txt`: Perl hash config wiring all VEP plugins above (paths, thresholds) into a single `--plugin` block

`resources_remapped/vep_gtf/vep_plugin_files.csv`: list of plugin Perl modules/config to stage as VEP `--dir_plugins` input, for `--vep_plugin_files`
