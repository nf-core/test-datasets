from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def reverse_quality(qual):
    return qual[::-1]

def generate_paired_end(input_fastq, output_fastq1, output_fastq2):
    with open(input_fastq, "r") as infile, open(output_fastq1, "w") as outfile1, open(output_fastq2, "w") as outfile2:
        for record in SeqIO.parse(infile, "fastq"):
            # Forward read
            SeqIO.write(record, outfile1, "fastq")
            
            # Reverse read (fake)
            reverse_seq = reverse_complement(str(record.seq))
            reverse_qual = reverse_quality(record.letter_annotations["phred_quality"])
            
            # Create a new SeqRecord for the reverse read
            reverse_record = SeqRecord(
                Seq(reverse_seq),
                id=record.id,
                description=record.description,
                letter_annotations={"phred_quality": reverse_qual}
            )
            
            SeqIO.write(reverse_record, outfile2, "fastq")



# Example usage
input_fastq = "Sample1_S1_L001_R1_001.fastq" #original single end fastq
output_fastq1 = "Sample1_R1.fastq" #name of output R1 fastq
output_fastq2 = "Sample1_R2.fastq" #name of output R2 fastq

generate_paired_end(input_fastq, output_fastq1, output_fastq2)

