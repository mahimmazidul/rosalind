# BA1N (d-Neighborhood)
# https://rosalind.info/problems/ba1n/
#
# What the algorithm is:
# Generate the d-neighborhood of a DNA string — all k-mers whose
# Hamming distance from the pattern does not exceed d.
#
# How it works conceptually:
# Recursively build all strings within distance d: at each position,
# either keep the original nucleotide (spending 0 mismatches) or
# substitute one of the other three (spending 1 mismatch), recursing
# on the remainder with the remaining mismatch budget.
#
# Time and space complexity:
# Time: O(neighborhood_size * k) — generate each neighbor string.
# Space: O(neighborhood_size * k) — store all neighbors.
#
# String representation: Recursive generation; output as a set of strings.
#
# Where this appears in bioinformatics:
# The d-neighborhood is the core search space of motif finding —
# finding a k-mer that appears (with up to d mismatches) in every
# string of a set of DNA sequences.

import sys
import os


def neighbors(pattern, d):
    """Generate all k-mers within Hamming distance d of pattern.

    Args:
        pattern: A DNA string.
        d: Maximum Hamming distance.

    Yields:
        All DNA strings of the same length within distance d of pattern.
    """
    if d == 0:
        yield pattern
        return
    if len(pattern) == 0:
        yield ''
        return

    first = pattern[0]
    suffix = pattern[1:]

    # Keep first character, spend mismatches on suffix
    for suffix_neighbor in neighbors(suffix, d):
        yield first + suffix_neighbor

    # Change first character (spend 1 mismatch), spend rest on suffix
    if d > 0:
        for base in 'ACGT':
            if base != first:
                for suffix_neighbor in neighbors(suffix, d - 1):
                    yield base + suffix_neighbor


def main():
    input_file = 'ba1n_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().strip().splitlines()
    else:
        lines = sys.stdin.read().strip().splitlines()

    pattern = lines[0].strip()
    d = int(lines[1].strip())

    result = list(neighbors(pattern, d))
    print(' '.join(result))


if __name__ == '__main__':
    main()
