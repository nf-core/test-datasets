Sample data from FASTP. This data was obtained running `nextflow run smrnaseq -profile test,docker -r dev --save_merged` with the following custom `publishDir`:

```
withName: '.*:FASTQ_FASTQC_UMITOOLS_FASTP:FASTP' {
        ext.args = [ "",
            params.trim_fastq                           ? "" : "--disable_adapter_trimming",
            params.clip_r1 > 0                          ? "--trim_front1 ${params.clip_r1}" : "", // Remove bp from the 5' end of read 1.
            params.three_prime_clip_r1 > 0              ? "--trim_tail1 ${params.three_prime_clip_r1}" : "", // Remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed.
            params.fastp_min_length > 0                 ? "-l ${params.fastp_min_length}" : "",
            params.fastp_max_length > 0                 ? "--max_len1 ${params.fastp_max_length}" : "",
            params.three_prime_adapter == "auto-detect" ?  "" : "--adapter_sequence ${params.three_prime_adapter}"
        ].join(" ").trim()
        publishDir = [
            [
                path: { "${params.outdir}/fastp/fastq" },
                mode: params.publish_dir_mode,
                pattern: "*.fastp.fastq.gz",
                enabled: params.save_merged
            ]
        ]
    }
```

Size was reduced with: `zcat Clone1_N3.fastp.fastq.gz | head -n 4000000 | gzip > small_Clone1_N3.fastp.fastq.gz`