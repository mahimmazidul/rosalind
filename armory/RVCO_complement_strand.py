# RVCO (Complementing a Strand of DNA)
# https://rosalind.info/problems/rvco/
#
# What a DNA complement is:
# DNA is double-stranded, with the two strands running in opposite
# directions (antiparallel). The bases on opposite strands always pair
# the same way: A with T, C with G. The complement of a base is the base
# it pairs with on the opposite strand.
#
# What reverse complement means:
# The reverse complement of a DNA string is formed by reversing the
# string and then taking the complement of each base. For example, the
# reverse complement of "GTCA" is "TGAC" (reverse: "ACTG", complement:
# "TGAC"). This matters because either strand of DNA can be the coding
# strand, so a sequence and its reverse complement represent the same
# double-stranded DNA region.
#
# Why this matters in real biology:
# When analyzing genomic DNA, researchers must consider both strands.
# Restriction enzymes, promoters, and other functional elements can
# appear on either strand. A sequence that equals its own reverse
# complement is a reverse palindrome — these are important because many
# restriction enzymes recognize and cut at such palindromic sites.
#
# What this solution does:
# It reads a collection of DNA strings in FASTA format, computes the
# reverse complement of each using Biopython's Seq object, and counts
# how many strings equal their own reverse complement.


from Bio.Seq import Seq
from Bio import SeqIO
import sys
import os


def count_reverse_palindromes(input_path):
    """Count DNA strings that equal their reverse complement.

    Args:
        input_path: Path to a FASTA file containing DNA strings.

    Returns:
        The number of strings that match their reverse complement.
    """
    records = SeqIO.parse(input_path, 'fasta')
    count = 0
    for record in records:
        seq = str(record.seq)
        rc = str(Seq(seq).reverse_complement())
        if seq == rc:
            count += 1
    return count


def main():
    input_file = 'rvco_input.txt'
    if os.path.exists(input_file):
        result = count_reverse_palindromes(input_file)
    else:
        import tempfile
        text = sys.stdin.read()
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(text)
            temp_path = f.name
        try:
            result = count_reverse_palindromes(temp_path)
        finally:
            os.unlink(temp_path)

    print(result)


if __name__ == '__main__':
    main()
