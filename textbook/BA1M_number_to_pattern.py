# BA1M (NumberToPattern)
# https://rosalind.info/problems/ba1m/
#
# What the algorithm is:
# Convert an integer index to its corresponding DNA string of length k,
# using base-4 encoding (0=A, 1=C, 2=G, 3=T).
#
# How it works conceptually:
# Repeatedly divide index by 4, collecting remainders as base-4 digits
# (least significant first). Reverse and map digits to nucleotides.
# Pad with leading A's (digit 0) to reach length k.
#
# Time and space complexity:
# Time: O(k) — k divisions by 4.
# Space: O(k) — output string.
#
# String representation: Base-4 digit mapping (inverse of PatternToNumber).
#
# Where this appears in bioinformatics:
# NumberToPattern is used when iterating over all possible k-mers in
# lexicographic order — e.g., enumerating candidate motifs in motif
# finding algorithms.

import sys
import os


def number_to_pattern(index, k):
    """Convert an integer index to a DNA string of length k.

    Args:
        index: An integer in [0, 4^k - 1].
        k: Length of the output DNA string.

    Returns:
        A DNA string of length k.
    """
    num_to_symbol = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    digits = []
    for _ in range(k):
        digits.append(num_to_symbol[index % 4])
        index //= 4
    return ''.join(reversed(digits))


def main():
    input_file = 'ba1m_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().strip().splitlines()
    else:
        lines = sys.stdin.read().strip().splitlines()

    index = int(lines[0].strip())
    k = int(lines[1].strip())

    print(number_to_pattern(index, k))


if __name__ == '__main__':
    main()
