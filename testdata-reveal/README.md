# Test datasets

* `sc5p_v2_hs_PBMC_1k_b_airr_rearrangement.tsv`

VDJ Ig - AIRR rearrangement (TSV) from Single Cell Immune Profiling Dataset by 
Cell Ranger 4.0.0, Human PBMC from a Healthy Donor, 1k cells (v2)

URL: https://support.10xgenomics.com/single-cell-vdj/datasets/4.0.0/sc5p_v2_hs_PBMC_1k

Published on August 25, 2020. This dataset is licensed under the Creative 
Commons Attribution license.

* `sc5p_v2_mm_c57bl6_splenocyte_1k_b_airr_rearrangement.tsv`

VDJ Ig - AIRR rearrangement (TSV) from Single Cell Immune Profiling Dataset by 
Cell Ranger 4.0.0,Splenocytes from C57BL/6 mice, 1k cells (v2)

URL: https://support.10xgenomics.com/single-cell-vdj/datasets/4.0.0/sc5p_v2_mm_c57bl6_splenocyte_1k

Published on August 25, 2020. This dataset is licensed under the Creative Commons Attribution license.

* `sc5p_v2_hs_B_prevax_10k_5gex_B_vdj_b_airr_rearrangement.tsv`
VDJ Ig - AIRR rearrangement (TSV). Human B cells isolated from a healthy female donor 
pre-influenza vaccination were obtained by 10x Genomics from Discovery Life Sciences.

URL: https://support.10xgenomics.com/single-cell-vdj/datasets/5.0.0/sc5p_v2_hs_B_prevax_10k_5gex_B

Single Cell Immune Profiling Dataset by Cell Ranger 5.0.0
Published on January 7, 2021
This dataset is licensed under the Creative Commons Attribution license.

* `sc5p_v2_hs_B_postvax_10k_5gex_B_vdj_b_airr_rearrangement.tsv`

VDJ Ig - AIRR rearrangement (TSV). Human B cells isolated from a healthy 
female donor 14 days post-influenza vaccination 
were obtained by 10x Genomics from Discovery Life Sciences.

URL: https://support.10xgenomics.com/single-cell-vdj/datasets/5.0.0/sc5p_v2_hs_B_postvax_10k_5gex_B

Single Cell Immune Profiling Dataset by Cell Ranger 5.0.0
Published on January 7, 2021
This dataset is licensed under the Creative Commons Attribution license.

* `vquest_airr_sc5p_v2_hs_PBMC_1k_b_airr_rearrangement.tsv`

Sequences from `sc5p_v2_hs_PBMC_1k_b_airr_rearrangement.tsv` were stored in a fasta file and
submitted to IMGT/HighV-QUEST (2021-05-25), requesting results in the AIRR format. 

* `sc5p_v2_hs_PBMC_1k_b_airr_rearrangement_sequences_igblast.tsv`

Sequences from `sc5p_v2_hs_PBMC_1k_b_airr_rearrangement.tsv` were stored in a fasta file and
run through IgBLAST with with Immcantation 4.1.0, requesting results in the AIRR format:

```
docker run -v $(pwd):/data:z --workdir /data immcantation/suite:4.1.0 \
   AssignGenes.py igblast -s sc5p_v2_hs_PBMC_1k_b_airr_rearrangement_sequences.fasta \
   -b /usr/local/share/igblast \
   --organism human --loci ig --format airr --nproc 1
```

* `bulk-Laserson-2014.fasta`

5,000 processed reads from a bulk sequencing experiment. Reads belong to one healthy donor (PGP1) 3 weeks after flu vaccination (Laserson et al. (2014), DOI: 10.1073/pnas.1323862111). As part of the processing, each sequence has been annotated with the isotype (`C_CALL`).
