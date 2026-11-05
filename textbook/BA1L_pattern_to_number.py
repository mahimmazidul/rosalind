# BA1L (PatternToNumber)
# https://rosalind.info/problems/ba1l/
#
# What the algorithm is:
# Convert a DNA string to a number by treating it as a base-4 integer
# with A=0, C=1, G=2, T=3.
#
# How it works conceptually:
# Each nucleotide is a base-4 digit. Reading left to right, the number
# is accumulated as number = number * 4 + digit_value.
#
# Time and space complexity:
# Time: O(k) — one pass over the string.
# Space: O(1) — just the accumulator.
#
# String representation: Base-4 digit mapping.
#
# Where this appears in bioinformatics:
# PatternToNumber is the core of k-mer indexing — mapping k-mers to
# array positions for frequency arrays, k-mer spectra, and hash-based
# sequence lookup.

import sys
import os


def pattern_to_number(pattern):
    """Convert a DNA k-mer to its base-4 index.

    Args:
        pattern: A DNA string over {A, C, G, T}.

    Returns:
        An integer in [0, 4^k - 1].
    """
    symbol_to_num = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    number = 0
    for base in pattern:
        number = number * 4 + symbol_to_num[base]
    return number


def main():
    input_file = 'ba1l_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().strip().splitlines()
    else:
        lines = sys.stdin.read().strip().splitlines()

    pattern = lines[0].strip()
    print(pattern_to_number(pattern))


if __name__ == '__main__':
    main()
