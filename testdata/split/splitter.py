import gzip

input_file = "C1-N1-R1_S4_L001_R1_001.fastq.gz"
output_file1 = "C1-N1-R1_S4_L001_R1_001_split1.fastq.gz"
output_file2 = "C1-N1-R1_S4_L001_R1_001_split2.fastq.gz"

with gzip.open(input_file, "rt") as infile, \
     gzip.open(output_file1, "wt") as outfile1, \
     gzip.open(output_file2, "wt") as outfile2:

    record = []
    for i, line in enumerate(infile):
        record.append(line)
        if len(record) == 4:  # A complete FASTQ record
            if (i // 4) % 2 == 0:
                outfile1.writelines(record)
            else:
                outfile2.writelines(record)
            record = []  # Reset for the next record
