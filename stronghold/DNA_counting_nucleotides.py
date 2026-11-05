# DNA (Counting DNA Nucleotides)
# https://rosalind.info/problems/dna/
#
# What is DNA?
# DNA is the molecule that carries the genetic instructions for the
# development, growth, and reproduction of all known living organisms. It is
# made up of two long strands that twist around each other in a double helix.
# Each strand is a chain of smaller building blocks called nucleotides.
#
# What are nucleotides?
# A nucleotide is one of the basic units of DNA. There are four types, each
# identified by a single letter:
#   A - Adenine
#   C - Cytosine
#   G - Guanine
#   T - Thymine
# The order of these letters along a DNA strand encodes genetic information,
# much like letters in a sentence.
#
# Why counting nucleotides matters in real biology:
# - Sequencing quality: When scientists sequence DNA, they check the counts of
#   each nucleotide to spot errors or biases in the sequencing run.
# - GC content: The proportion of G and C nucleotides in a sequence (GC
#   content) is a quick measure of DNA stability, since G-C pairs are held
#   together by three hydrogen bonds while A-T pairs have only two.
# - Mutation studies: Comparing nucleotide counts between a healthy sequence
#   and a mutated one can reveal insertions, deletions, or substitutions.
#
# This is Rosalind problem #1 in the Stronghold track. The task is simple:
# given a DNA string, count how many times A, C, G, and T appear, and output
# the four counts separated by spaces in that order.


def count_nucleotides(dna):
    """Count the occurrences of A, C, G, and T in a DNA string.

    Args:
        dna: A string consisting of the characters 'A', 'C', 'G', and 'T'.

    Returns:
        A dictionary with keys 'A', 'C', 'G', 'T' and integer counts.
    """
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for base in dna:
        if base in counts:
            counts[base] += 1
    return counts


def main():
    import sys
    import os

    input_file = 'input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            dna = f.read().strip()
    else:
        dna = sys.stdin.read().strip()

    counts = count_nucleotides(dna)
    print(f"{counts['A']} {counts['C']} {counts['G']} {counts['T']}")


if __name__ == '__main__':
    main()
