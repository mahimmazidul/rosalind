# BA1C (Reverse Complement)
# https://rosalind.info/problems/ba1c/
#
# What the algorithm is:
# Reverse complement of a DNA string — reverse the string, then replace
# each nucleotide with its complement (A↔T, C↔G).
#
# How it works conceptually:
# DNA strands are antiparallel. The reverse complement represents the
# opposite strand read in the 5'→3' direction. Reverse the input string
# and map each base to its complement.
#
# Time and space complexity:
# Time: O(n) — one pass to reverse, one pass to complement.
# Space: O(n) — output string.
#
# String representation: Python string; translation via a lookup dict.
#
# Where this appears in bioinformatics:
# Reverse complement is used everywhere — designing PCR primers,
# searching both strands of DNA for a motif, reading genes on the
# minus strand, and comparing sequences regardless of orientation.

import sys
import os


def reverse_complement(pattern):
    """Return the reverse complement of a DNA string.

    Args:
        pattern: A DNA string over {A, C, G, T}.

    Returns:
        The reverse complement string.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[b] for b in reversed(pattern))


def main():
    input_file = 'ba1c_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().strip().splitlines()
    else:
        lines = sys.stdin.read().strip().splitlines()

    pattern = lines[0].strip()
    print(reverse_complement(pattern))


if __name__ == '__main__':
    main()
