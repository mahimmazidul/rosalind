# TFSQ (FASTQ Format Introduction)
# https://rosalind.info/problems/tfsq/
#
# What FASTQ is:
# FASTQ is the standard format for storing high-throughput sequencing
# output. It extends FASTA by adding per-base quality scores. Each record
# has four lines: a header line starting with @, the nucleotide sequence,
# a line starting with + (often repeating the header or empty), and a line
# of quality scores encoded as ASCII characters.
#
# What Phred quality scores are:
# Sequencing machines make errors. The Phred quality score Q quantifies
# the confidence in each base call: Q = -10 * log10(P), where P is the
# probability that the base is called incorrectly. Q=10 means 90% accuracy,
# Q=20 means 99%, Q=30 means 99.9%. In FASTQ files, the Phred score is
# encoded as an ASCII character by adding an offset of 33 (so Q=10 becomes
# '+', ASCII 43). Higher scores mean higher confidence.
#
# Why quality scores matter:
# Quality scores let researchers filter out low-quality reads, trim
# poor-quality ends, and weight bases appropriately in downstream analysis
# like variant calling or assembly. Without quality information, a
# sequencing error looks identical to a real genetic variant.
#
# What this solution does:
# It reads a FASTQ file using Biopython's SeqIO, which handles the format
# parsing automatically, and writes each record as FASTA (header with >
# instead of @, sequence only, dropping quality scores).


from Bio import SeqIO
import sys
import os


def fastq_to_fasta(input_path, output_path=None):
    """Convert a FASTQ file to FASTA format.

    Reads FASTQ records from input_path and writes them as FASTA to
    output_path (or stdout if output_path is None).

    Args:
        input_path: Path to the input FASTQ file.
        output_path: Optional path to the output FASTA file.
    """
    records = SeqIO.parse(input_path, 'fastq')
    if output_path:
        SeqIO.write(records, output_path, 'fasta')
    else:
        for record in records:
            print(f'>{record.id}')
            print(record.seq)


def main():
    input_file = 'tfsq_input.txt'
    if os.path.exists(input_file):
        fastq_to_fasta(input_file)
    else:
        # Read FASTQ from stdin
        import tempfile
        text = sys.stdin.read()
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
            f.write(text)
            temp_path = f.name
        try:
            fastq_to_fasta(temp_path)
        finally:
            os.unlink(temp_path)


if __name__ == '__main__':
    main()
