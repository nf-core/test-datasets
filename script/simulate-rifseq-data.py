#!/usr/bin/env python

import argparse
import random
from itertools import groupby

# Argument parsing
parser = argparse.ArgumentParser(description = 'Simulate a RIF-Seq dataset.')
parser.add_argument('fasta',
                    type    = str,
                    metavar = 'fasta',
                    help    = 'A reference FASTA to create reads from')
parser.add_argument('output',
                    type    = str,
                    metavar = 'output',
                    help    = 'Output path for the resulting FASTQ file')
parser.add_argument('-l', '--length',
                    type    = int,
                    dest    = "length",
                    default = 76,
                    metavar = '',
                    help    = 'The read length of to simulate')
parser.add_argument('-n', '--number',
                    type    = int,
                    dest    = "number",
                    default = 1990,
                    metavar = '',
                    help    = 'Number of reads to generate per barcode')
parser.add_argument('-u', '--unknown_number',
                    type    = int,
                    dest    = "unknown_number",
                    default = 10,
                    metavar = '',
                    help    = 'Number of unknown reads to generate per barcode')
parser.add_argument('-s', '--seed',
                    type    = int,
                    dest    = "seed",
                    default = 42,
                    metavar = '',
                    help    = 'The seed to use for random number generation')
args = parser.parse_args()


def extract_fasta_information(fasta_file):
    """ Extract sequenced and IDs from an Ensembl FASTA file

    Extracts both FASTA sequences and their transcript IDs from an Ensembl
    transcript FASTA file, returning a dict with the former as values and the
    latter as keys.

    Modified from Brent Pedersen
    "Correct Way To Parse A Fasta File In Python"
    https://www.biostars.org/p/710/
    """

    # Initialise dict
    fasta_dict = dict()

    # Open the file
    fh = open(fasta_file)

    # Ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate
    fasta_iter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    # Iterate over each line
    for line in fasta_iter:

        # Get the transcript ID
        transcript = next(line)[1:].strip().split()[0]

        # Join all sequence lines to one
        sequence = ''.join(seq.strip() for seq in next(fasta_iter))

        # Add to dict
        fasta_dict[transcript] = sequence

    # Return final FASTA dictionary
    return(fasta_dict)


def check_transcript_length(fasta_dict, length):
    """ Select and check a random FASTA record for sequence length

    Selects a random record from a FASTA dictionary and checks the length of the
    corresponding sequence, to see if it can be used to create a read with at
    least the specified length.
    """

    # Get transcript ID and sequence
    transcript = random.choice(list(fasta_dict.keys()))
    sequence = fasta_dict[transcript]

    # Check if transcript length is sufficient for building a read
    if len(sequence) < length - 18:
        return None
    else:
        return [transcript, sequence]


def build_record(transcript, sequence, barcode, seq_length):
    """ Build a RIF-Seq FASTQ record

    Build a complete RIF-Seq FASTQ record from input transcript ID, sequence,
    RIF-Seq barcode and specified sequence length.
    """

    # Build read sequence from barcode, random decamer and sequence
    start_pos = random.randrange(0, len(sequence) - (seq_length - 18))
    random_decamer = ''.join(random.choices("ATCG", k = 10))
    read = barcode + \
        random_decamer + \
        sequence[start_pos:start_pos + seq_length - 18]

    # Build random quality scores
    quality_scores = ''.join(random.choices(qual, k = seq_length))

    # Build header
    header = '@' + transcript + \
            ' barcode:' + barcode + \
            ' decamer:' + random_decamer

    # Build entire FASTQ record and return
    record = header + '\n' + read + '\n+\n' + quality_scores + '\n'
    return(record)


# Set seed
#  random.seed(42)
random.seed(args.seed)

# Define possible values for FASTQ quality scores
qual = '!#$%&()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ' + \
    '[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'

# Load FASTA
fasta_dict = extract_fasta_information(args.fasta)

# Define barcodes: 4 genuine RIF-Seq sequences and 1 non-RIF-Seq sequence
barcode_sequences = ['TCATGAAT', 'GAACACTT', 'CAAGTAAC', 'ACCGACAC', 'AATTCCGG']

# Open output file for writing
fh = open(args.output, "w")

# Create `number` of FASTQ recrods for each barcode
for barcode in barcode_sequences:

    # Progress
    print('Building records for RIF-Seq barcode ' + barcode + ' ...')
    nn = 0

    # Define number of reads to get depending on barcode sequence
    if barcode == "AATTCCGG":
        num_reads = args.unknown_number * 4
    else:
        num_reads = args.number

    # Build `number` of RIF-Seq records
    while nn < num_reads:

        # Get transcript ID/sequence and build FASTQ record
        info = check_transcript_length(fasta_dict, args.length)
        if info != None:
            record = build_record(info[0], info[1], barcode, args.length)
            fh.write(record)
            nn += 1
        else:
            continue
