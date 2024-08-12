# SingleCell-RNA_P3_2 Dataset

## Overview

This repository contains a processed dataset "SingleCell-RNA_P3_2," obtained from Illumina BaseSpace. The dataset includes a flowcell tar.gz file and the SampleSheet.csv file included in this folder.

## Dataset Processing

1. **Download**: The dataset was downloaded using [`basemount`](https://basemount.basespace.illumina.com/#install).
2. **Trimming**: cbcl files were trimmed to include only the first tile using [cbcltrimmer/trim.py](https://github.com/Aratz/cbcltrimmer/blob/main/trim.py).
3. **RunInfo.xml Update**: The `TileCount` was changed to `TileCount="1"`.
4. **Tile Entries Simplification**: Only the first tile entries (1_1101 and 2_1101) were retained.
