# Pixelator 0.29 test data overlay

Updated fixtures for pixelator 0.29.0 live under `new-test-data/pna/`.
Unchanged files are symlinks to `testdata/pna/`; changed or new files are stored here.

The v2 samplesheet with the prerelease panel is at
`new-test-data/samplesheet/pna/samplesheet_proxiome_v2.csv`.

When running nf-core/pixelator tests locally, point to this repo clone so symlinks resolve:

```bash
export NFT_TESTDATA_BASE_PATH="$PWD/"
```
