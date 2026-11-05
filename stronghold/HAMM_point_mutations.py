# HAMM (Counting Point Mutations)
# https://rosalind.info/problems/hamm/
#
# What a point mutation is:
# A point mutation is a change in a single nucleotide of a DNA sequence. It
# is the smallest possible change to the genetic code — one letter swapped
# out for another at a specific position in the genome.
#
# The three types of point mutations:
# Substitution: one base is replaced by another (e.g., A becomes G). This is
# what Hamming distance measures — it counts how many positions differ
# between two equal-length strings. Insertion: an extra base is added,
# shifting everything after it. Deletion: a base is removed, again shifting
# everything after it. Insertions and deletions change the length of the
# sequence, so they can't be measured by Hamming distance, which requires
# both strings to be the same length.
#
# Why point mutations matter:
# Not all point mutations have the same effect. Some are silent — the codon
# still codes for the same amino acid because the genetic code is redundant.
# Some are missense — the codon now codes for a different amino acid, which
# can change the protein's shape and function. Some are nonsense — the
# codon becomes a stop codon, truncating the protein entirely. A famous
# example: sickle cell anaemia is caused by a single point mutation where
# adenine is replaced by thymine in the beta-globin gene, changing one amino
# acid in hemoglobin and distorting red blood cells into a sickle shape.
#
# Hamming distance as a concept:
# Richard Hamming developed this idea in information theory as a way to
# measure how many positions differ between two equal-length strings. It was
# originally used for error detection and correction in telecommunications
# — if a received message differs from a valid codeword by a Hamming distance
# of 1, you know a single-bit error occurred. Molecular biologists adopted
# it to compare DNA sequences because point mutations are exactly the kind of
# single-position errors that Hamming distance captures.
#
# Where this is used in practice:
# - Comparing sequencing reads to a reference genome to find variants.
# - Measuring genetic divergence between species or individuals.
# - Error correction in high-throughput sequencing, where reads are compared
#   to known sequences to identify and fix sequencing mistakes.


def hamming_distance(s, t):
    """Calculate the Hamming distance between two equal-length strings.

    The Hamming distance is the number of positions at which the
    corresponding symbols differ. It only applies to strings of the same
    length, since it measures point substitutions.

    Args:
        s: First DNA string.
        t: Second DNA string (same length as s).

    Returns:
        The number of differing positions between s and t.
    """
    if len(s) != len(t):
        raise ValueError("Strings must be of equal length for Hamming distance")
    return sum(1 for a, b in zip(s, t) if a != b)


def main():
    import sys
    import os

    input_file = 'hamm_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().splitlines()
    else:
        lines = sys.stdin.read().splitlines()

    # Filter out empty lines and take the first two non-empty lines
    dna1 = None
    dna2 = None
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if dna1 is None:
            dna1 = line
        elif dna2 is None:
            dna2 = line
            break

    print(hamming_distance(dna1, dna2))


if __name__ == '__main__':
    main()
