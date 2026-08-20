# ONT pod5 test data — Homo sapiens

## HG002_PAW70337_giab_10reads.pod5

| Field     | Value                                                                                              |
| --------- | -------------------------------------------------------------------------------------------------- |
| Sample    | HG002 (Genome in a Bottle reference sample)                                                        |
| Run       | PAW70337                                                                                           |
| Chemistry | 5kHz, R10.4.1                                                                                      |
| Reads     | 10 (subset)                                                                                        |
| Size      | ~2.2 MB                                                                                            |
| Source    | `s3://ont-open-data/giab_2025.01/flowcells/HG002/PAW70337/pod5/PAW70337_66b2eea5_de8117b1_33.pod5` |
| Access    | Public, no authentication required (`--no-sign-request`)                                           |
| Released  | 2025-01-14                                                                                         |

### How this file was created

```bash
# Download original (530 MB)
aws s3 cp --no-sign-request \
  s3://ont-open-data/giab_2025.01/flowcells/HG002/PAW70337/pod5/PAW70337_66b2eea5_de8117b1_33.pod5 .

# Extract 10 reads
pod5 view PAW70337_66b2eea5_de8117b1_33.pod5 --no-header | head -10 | \
  awk '{print "HG002_PAW70337_giab_10reads.pod5," $1}' > subset.csv
echo "target_filename,read_id" | cat - subset.csv > subset_with_header.csv
pod5 subset PAW70337_66b2eea5_de8117b1_33.pod5 --csv subset_with_header.csv --output .
```

### Usage in nf-core modules tests

```groovy
input[0] = [
    [ id: 'HG002' ],
    file("${params.modules_testdata_base_path}genomics/homo_sapiens/nanopore/pod5/HG002_PAW70337_giab_10reads.pod5",
         checkIfExists: true)
]
```
