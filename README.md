# test-datasets: `pixelator`

This branch contains test data to be used for automated testing with the [nf-core/pixelator](https://github.com/nf-core/pixelator) pipeline.

## Content of this repository

- `testdata/` :
    - `mpx/modules` : Testdata for pipeline local module tests
    - `mpx/fastq` : Input fastq files for pipeline level tests
- `panels/`: panel files to test passing custom panels

- `samplesheet/samplesheet.csv`: Experiment design file for minimal test dataset.
- `samplesheet/samplesheet_full.csv`: Experiment design file for full test dataset.
- `samplesheet/samplesheet_v2.csv`: Experiment design file for minimal test dataset using the v2 panel.
- `samplesheet/samplesheet_mpx_scsp_v1_immunology1.csv` : Experiment design file for minimal test dataset using the same data as the local module tests

## Test datasets origin

Molecular pixelation data retrieved from public datasets provided by Pixelgen Technologies AB.

- https://software.pixelgen.com/datasets/1k-human-pbmcs-v1.0-immunology-I
- https://software.pixelgen.com/datasets/uropod-t-cells-v1.0-immunology-I

Immunology II datasets are taken from:

s3://pixelgen-technologies-datasets/mpx-datasets/scsp/2.0/5-donors-pbmcs-v2.0/


### Sampling procedure


Data was subsampled using `seqkit sample` to 200k or 300k paired-end reads.
The default random seed was used.

```bash
seqkit sample -2 -n 300000 Uropod_control_R1_001.fastq.gz -o uropod_300k_control_R1_001.fastq.gz
seqkit sample -2 -n 300000 Uropod_control_R2_001.fastq.gz -o uropod_300k_control_R2_001.fastq.gz
seqkit sample -2 -n 200000 Sample01_human_pbmcs_unstimulated_200k_R1_001.fastq.gz -o Sample01_human_pbmcs_unstimulated_200k_R1_001.fastq.gz
seqkit sample -2 -n 200000 Sample01_human_pbmcs_unstimulated_200k_R2_001.fastq.gz -o Sample01_human_pbmcs_unstimulated_200k_R2_001.fastq.gz
seqkit sample -2 -n 500000 Sample01_PBMC_1085_r1_R1_001.fastq.gz -o v2/Sample01_PBMC_1085_500k_r1_R1_001.fastq.gz
seqkit sample -2 -n 500000 Sample01_PBMC_1085_r1_R2_001.fastq.gz -o v2/Sample01_PBMC_1085_500k_r1_R2_001.fastq.gz
```

To test input sample concatenation the uropod_control sample was split using `seqkit split2`

```bash
seqkit split2 -p 2 --read1 uropod_control_300k_R1_001.fastq.gz --read2 uropod_control_300k_R2_001.fastq.gz
```

Seqkit is available on bioconda.


## `testdata/mpx`

Data in `testdata/mpx` has been generated for nf-test local module and pipeline tests.

The raw data and processed data are downloaded from the ["1k human PBMCs v1 dataset"](https://software.pixelgen.com/datasets/1k-human-pbmcs-v1.0-immunology-I)

- Sample01_human_pbmcs_unstimulated.layout.dataset.pxl
- Sample01_human_pbmcs_unstimulated_R1_001.fastq.gz
- Sample01_human_pbmcs_unstimulated_R2_001.fastq.gz


### Subsampling procedure

Data is subsampled by picking a few cells from a processed pixeldataset file and creating a list of umi + upi-a sequences for all reads in those cells.
Next we collect the sequence IDs in the raw input data for presence of these UMI+UPI-A combinations in the respective positions.
This procedure guarantees that the subsampled dataset will still generate a few well behaved cells.

The raw data is then filtered with the fastq ids and further subsampled to 100k reads using seqkit.

```python
import numpy as np
import dnaio
from pixelator import read
import itertools
import polars as pl
import pandas as pd
from pathlib import Path
from pixelator.utils import reverse_complement
from pixelator.config import config as pixelator_config
from pixelator.config.assay import get_position_in_parent


def filter_records(source, target, whitelist, rev_comp=False):
    records_handled = 0
    records_matched = 0

    assay = pixelator_config.get_assay("D21")
    umib_start, umib_end = get_position_in_parent(assay, "umi-b")
    upia_start, upia_end = get_position_in_parent(assay, "upi-a")

    amplicon = assay.get_region_by_id("amplicon")

    r2_umib_start, r2_umib_end = amplicon.min_len - umib_end, amplicon.min_len - umib_start
    r2_upia_start, r2_upia_end = amplicon.min_len - upia_end, amplicon.min_len - upia_start

    print(r2_umib_start, r2_umib_end)

    with open(target, "w") as out:
        with dnaio.open(source, mode='r') as f:
            for record in f:
                umib = record.sequence[r2_umib_start:r2_umib_end]
                upia = record.sequence[r2_upia_start:r2_upia_end]

                if rev_comp:
                    umib = reverse_complement(umib)
                    upia = reverse_complement(upia)

                seq_match = umib + upia
                if seq_match in whitelist:
                    records_matched += 1
                    out.write(record.id + '\n')
            
                records_handled += 1

                if records_handled % 1_000_000 == 0:
                    print(f"Records handled: {records_handled}, matched: {records_matched}")

# Read input file
pxl_file_path = Path.cwd() / "Sample01_human_pbmcs_unstimulated.layout.dataset.pxl"
R1 =  Path.cwd() / "Sample01_human_pbmcs_unstimulated_R1_001.fastq.gz"
R2 =  Path.cwd() / "Sample01_human_pbmcs_unstimulated_R2_001.fastq.gz"

# Partition by component
pxl_file = read(pxl_file_path)
df = pxl_file.edgelist_lazy.collect()
parts = df.partition_by("component", as_dict=True)

# Sort on component size
sorted_parts = np.array(sorted(parts.items(), key=(lambda t: t[1].shape[0]), reverse=True), dtype=object)

# Subset last 10 components and join the data frame
indices = np.random.randint(0, len(sorted_parts), 5)
component_subset = pl.concat(v for k, v in sorted_parts[indices])

# Collect umis
unique_umi_upia = component_subset[["umi", "upia"]].unique()
unique_seqs = set(unique_umi_upia.to_pandas().apply(lambda x: x["umi"] + x["upia"], axis=1))

# Create file with filtered ids
filtered_ids = Path.cwd() / "filtered_ids.txt"
filter_records(R2, filtered_ids, unique_seqs, rev_comp=True)
```

```bash
seqkit grep -f filtered_ids.txt Sample01_human_pbmcs_unstimulated_R1_001.fastq.gz -o Sample01_human_pbmcs_unstimulated_R1_001.subset.fastq.gz
seqkit sample -n 100000 Sample01_human_pbmcs_unstimulated_R1_001.subset.fastq.gz -o sample01_1k_pbmcs_scsp_v1_immunology1_R1.fastq.gz

seqkit grep -f filtered_ids.txt Sample01_human_pbmcs_unstimulated_R2_001.fastq.gz -o Sample01_human_pbmcs_unstimulated_R2_001.subset.fastq.gz
seqkit sample -n 100000 Sample01_human_pbmcs_unstimulated_R2_001.subset.fastq.gz -o sample01_1k_pbmcs_scsp_v1_immunology1_R2.fastq.gz
```