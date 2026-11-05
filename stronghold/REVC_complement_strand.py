# REVC (Complementing a Strand of DNA)
# https://rosalind.info/problems/revc/
#
# What the two strands of DNA are, and why they run antiparallel:
# DNA is a double helix made of two strands that wrap around each other.
# Each strand has a direction: one runs 5' to 3' (five-prime to three-prime)
# and the other runs 3' to 5'. The numbers refer to the carbon atoms on the
# sugar molecule in the backbone where the phosphate group attaches. Because
# the two strands are oriented in opposite directions, they are antiparallel.
# This matters because DNA polymerase and other enzymes can only read a strand
# in one direction — 5' to 3' — so the complementary strand must be read
# backwards relative to the original to make sense.
#
# Watson-Crick base pairing and why these specific pairs form:
# The two strands are held together by hydrogen bonds between nitrogenous
# bases. Adenine (A) always pairs with Thymine (T), and Guanine (G) always
# pairs with Cytosine (C). A and T form two hydrogen bonds between them,
# while G and C form three. That extra hydrogen bond makes G-C pairs
# stronger — it takes more energy to pull them apart — which is why
# GC-rich regions of DNA are more thermally stable and melt at higher
# temperatures than AT-rich regions.
#
# What reverse complement means biologically:
# If you have a DNA sequence written 5' to 3', its complementary strand runs
# 3' to 5' underneath it. To read that complementary strand in the
# biologically correct 5' to 3' direction, you have to read it backwards
# while swapping each base for its complement (A↔T, C↔G). That backward,
# base-swapped sequence is the reverse complement. For example, if the
# original strand is AAAACCCGGT, the complement is TTTTGGCCA, and read
# backwards (5' to 3') it becomes ACCGGGTTTT.
#
# Where this actually matters:
# PCR primer design relies on reverse complements — primers are short DNA
# sequences that bind to the complementary strand so polymerase can copy the
# target region. Restriction enzymes often recognize palindromic sequences,
# which read the same on both strands when each is read 5' to 3'. Finding
# these sites means computing reverse complements. So this simple-looking
# operation shows up constantly in molecular biology lab work.


def reverse_complement(dna):
    """Return the reverse complement of a DNA string.

    Each base is replaced by its Watson-Crick complement (A<->T, C<->G),
    and the resulting string is reversed so it reads 5' to 3'.

    Args:
        dna: A string of 'A', 'C', 'G', 'T' characters.

    Returns:
        The reverse complement string.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(dna))


def main():
    import sys
    import os

    input_file = 'revc_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            dna = f.read().strip()
    else:
        dna = sys.stdin.read().strip()

    print(reverse_complement(dna))


if __name__ == '__main__':
    main()
