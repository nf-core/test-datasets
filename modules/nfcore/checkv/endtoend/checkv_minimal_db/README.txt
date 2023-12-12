## Verion 1.5 (Jan 10, 2023)

The main change in v1.5 was to remove putatively non-viral sequences to minimize high-confidence matches to the CheckV database for non-viral sequences.

first, viral features and taxonomy were obtained for all sequences
* ran geNomad v1.3.1, viralVerify v1.1, CheckV v1.1
* see "genome_db/checkv_info.tsv" for details

a "viral score" was calculated by taking the sum of the following:
* checkv host genes > viral genes or checkv viral genes == 0: -1
* genomad virus score < 0.7: -1
* genomad 0 viral hallmarks: -1
* genomad >=1 plasmid hallmarks: -1
* viralverify Plasmid label: -1
* viralverify Chromosome label: -1
* checkv viral genes > 3 x host genes and viral genes >= 2: +1
* genomad score > 0.95: +1
* genomad viral hallmarks >=2: +1
* viralverify Virus label: +1

putatively non-viral sequences were removed as follows:
* having a score of -2 or less
* flagged as non-viral by at least 2/3 tools
* not classified as viral by any tool
* annotated as Caudoviricetes but having a genome size <10 Kbp
* exceptions to the rules above were made for Inovirus, RNA viruses from the RVMT, and RefSeq/ICTV genomes
* overall 1913 viral sequences were removed from the v1.4 database
* see "genome_db/changelog.tsv" for details

## Verion 1.4 (Aug 27, 2022)
* excluded genomes from additional NCBI bioprojects: "PRJNA393166", "PRJNA492716", "PRJEB23154", "PRJNA579886", "PRJDB11069"

## Verion 1.3 (Aug 21, 2022)
added new metagenome-assembed complete genomes:
* Lavidaviridae (ITRs, n=782, 17-32kbp, >=3 geNomad MGs)
* Inoviridae (DTRs, n=673, 4-14kbp, derived from PMID 31332386); replaces the 988 in v1.1
* Large Caudovirices w/o matches to the CheckV database (DTRs, n=3448, >50kbp, >=3 geNomad MGs, >=1 terminase/MCP/portal gene)
* NCLDV (DTRs, n=80, >100kbp, >=10 geNomad MGs)
* Putative complete RNA viruses from RVMT dataset (https://doi.org/10.1101/2022.02.15.480533)
* Asgard archaeal viruses (n=2, MDVT01000002, MDVT01000005)
* Environmental DJR phages with DTRs/ITRs (https://doi.org/10.1186/s12985-018-0974-y)

added new NCBI complete genomes:
* ANI < 95 or AF < 85 versus existing CheckV genomes
* genome_rep = "Full" & assembly_level in ("Complete genome", "Chromosome")
* genome length > 1500 bp
* excluded genomes from specific NCBI bioprojects with low quality genomes
* excluded genomes classified as proviruses
* excluded genomes with high kmer freq (>1.2)

additionally:
* rebuilt protein database using gene calls made with prodigal-gv
* add new file indicating the data source of reference genomes: checkv-db-v1.3/genome_db/checkv_source.tsv
* removed a handful of DTR genomes from the initial CheckV database that appeared abnormally long

## Verion 1.2 (May 10, 2022)
* removed DIAMOND database to reduce space
* DIAMOND database is now built locally after download

## Version 1.1 (Mar 18, 2022)
* added 988 complete Inovirus genomes from Roux et al. (PMID 31332386)
* added 5541 complete phage genomes from INPHARED (https://github.com/RyanCook94/inphared)
* INPHARED genomes were filtered to be distinct from existing CheckV genomes
* both sets of genomes were screened by CheckV to remove likely genome fragments and non-viral sequences
* sequences files were updated with new seqs (checkv_reps.fna, checkv_reps.faa, checkv_reps.dmnd)
* metadata files were NOT updated with new seqs (checkv_circular.tsv, checkv_genbank.tsv, & checkv_clusters.tsv)
* new database release should improve completeness estimates for Inovirus and possibly other groups

## Version 1.0 (Feb 2, 2021)
* stable v1.0 release

## Version 0.6 (May 6, 2020)
* added hmm_db/genome_lengths.tsv
* publication release

genome_db
|
|-checkv_circular.tsv
  metadata for circular viral contigs
|-checkv_genbank.tsv
  metadata for genbank genomes
|-checkv_clusters.tsv
  information for non-redundant genomes
  these are clustered at 95% ANI over 85% the length of both genomes
  clustering performed using a greedy centroid-based algorithm
  the cluster representatives are used by CheckV for mapping
|-checkv_reps.faa
  proteins for representative genomes
|-checkv_reps.fna
  genomic sequences for representative genomes
  terminal repeat trimmed for circular viral contigs
|-checkv_reps.dmnd
  DIAMOND database of protein sequences
|-checkv_reps.tsv
  basic info on representative genomes
|-checkv_error.tsv
  lookup table used to report error and assign confidence levels to CheckV completeness estimates

hmm_db
|
|-checkv_hmms.tsv
  metadata for 15,959 CheckV HMMs
|-checkv_hmms
  directory of HMMs
  each file contains 200 HMMs
  HMMs split into multiple files for improved parallelization
|-genome_lengths.tsv
  genome length statistics for each HMM
