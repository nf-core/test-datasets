The fastq file in this repo is a subset of 10 reads extracted from raw read SRR19726169 (available from any INSDC database), representing a chicken gut metagenome sample sequenced with PacBio HiFi.
This is a good test file for any tool filtering adapters from PacBio HiFi reads, as 1 read out of 10 contains an adapter.

The original file was selected as some of its reads still contained an adapter, while some others didn't.
The first 9 reads have been cropped out of the original file with `seqkit`. They don't contain any adapters.
To identify a read with adapters, `HifiAdapterFilt` was launched on the original raw reads. The tenth read was selected among the ones filtered out by `HiFiAdapterFilt`.
