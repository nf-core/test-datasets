Test samplesheets for nf-core/pairgenomealign
=============================================

CI tests (`testsamplesheet.csv`)
--------------------------------

Minimal data tests for continuous integration (CI). The sample sheet contains
raw URLs to SARS-CoV-2 sequences from this repository's `modules` branch on
GitHub.

Full tests (`testsamplesheet_full.csv`)
---------------------------------------

Full, real-sized tests to stress the pipeline and demonstrate its practical
use. The sample sheet contains download URLs from the NCBI website for
mammalian genomes publicly available in the GenBank (GCA IDs) or RefSeq (GCF
IDs) databases. Note that depending on the user's location, the download may
fail when using Nextflow while succeeding with other tools such as `wget` or
`curl`. This is repository-side behavior, possibly as a consequence of attempts
to control costs amid the sharp increase in network traffic caused by AI
companies in 2026. At the moment, no failures have been observed when running
the full tests from the Amazon Web Services (AWS) cloud.

Small tests (`testsamplesheet_small.csv`)
-----------------------------------------

Small tests on unicellular eukaryotic genomes (NCBI/GenBank; see above) to
quickly demonstrate the use of the pipeline while producing more realistic
output than the CI tests.
