# Test KOfam database for nf-core/metatdenovo

A tiny `ko_list.gz` + `profiles.tar.gz` pair, in exactly the format `KOFAMSCAN_DOWNLOAD`
downloads from https://www.genome.jp/ftp/db/kofam/, used to exercise the pipeline's KofamScan
subworkflow in `nf-test` via the `--kofam_ko_list_url`/`--kofam_profiles_url` parameters, without
downloading the real ~1.5 GB `profiles.tar.gz`.

## How this was built (2026-07-26)

nf-core/metatdenovo's `-profile test` (megahit + prodigal) ORFs were run once against the real,
full KOfam database (profiles dated 2022-12-31) on a separate machine, giving a ground-truth
`besthit` table of every ORF's best-scoring KO. From that, the 20 KOs with the most significant
(`*`-flagged, i.e. score above the KO's own threshold) hits were picked — no arbitrary cutoff was
needed since the count distribution already tapers off naturally (16, 13, 5, 5, 4, 4, 4, 4, then a
long tail of 3s):

| KO | Count | Definition |
| --- | --- | --- |
| K21572 | 16 | starch-binding outer membrane protein, SusD/RagB family |
| K21573 | 13 | TonB-dependent starch-binding outer membrane protein SusC |
| K07497 | 5 | putative transposase |
| K07165 | 5 | transmembrane sensor |
| K23514 | 4 | RNA polymerase sigma-19 factor, ECF subfamily |
| K07483 | 4 | transposase |
| K03531 | 4 | cell division protein FtsZ |
| K01546 | 4 | potassium-transporting ATPase potassium-binding subunit |
| K24119 | 3 | NADH-dependent peroxiredoxin subunit C [EC:1.11.1.26] |
| K15633 | 3 | 2,3-bisphosphoglycerate-independent phosphoglycerate mutase [EC:5.4.2.12] |
| K09955 | 3 | uncharacterized protein |
| K08218 | 3 | MFS transporter, PAT family, beta-lactamase induction signal transducer AmpG |
| K07480 | 3 | insertion element IS1 protein InsB |
| K07407 | 3 | alpha-galactosidase [EC:3.2.1.22] |
| K07126 | 3 | uncharacterized protein |
| K07042 | 3 | probable rRNA maturation factor |
| K06142 | 3 | outer membrane protein |
| K04079 | 3 | molecular chaperone HtpG |
| K03590 | 3 | cell division protein FtsA |
| K02968 | 3 | small subunit ribosomal protein S20 |

The dominance of K21572/K21573 (a SusC/SusD starch-utilization pair) is consistent with
*Bacteroides fragilis* being a major community member in this test dataset — see
`test_data/diamond_taxonomy/README.md` for the taxonomic side of the same test assembly.

Unlike the DIAMOND taxonomy fixture, the HMM profiles themselves can't be reconstructed from the
pipeline's own ORFs — they're the real, unmodified KEGG-curated `.hmm` files for these 20 KOs,
downloaded directly by the pipeline author, plus the matching rows (header + these 20 KOs) from
the real `ko_list` file (needed for per-KO score thresholds). Re-running against these 20 profiles
reproduces real, correctly-thresholded hits (94 of 158 total hits score above threshold when
checked against the pipeline's own test ORFs).

## Regenerating or extending

Run `nf-core/metatdenovo`'s `tests/default.nf.test` to reproduce the ORFs, run the real KofamScan
database against them to get a fresh ground-truth `besthit` table, pick whichever KOs are needed,
and pull just those `profiles/K*.hmm` files plus the matching `ko_list` rows (header + one row per
chosen KO) from a real KOfam database download. Repackage as `gzip ko_list` and
`tar -czf profiles.tar.gz profiles/`.
