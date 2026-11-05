# BA1G (Hamming Distance)
# https://rosalind.info/problems/ba1g/
#
# What the algorithm is:
# Hamming distance — number of positions at which two equal-length
# strings differ.
#
# How it works conceptually:
# Compare the two strings character by character, counting mismatches.
#
# Time and space complexity:
# Time: O(n) — single pass over both strings.
# Space: O(1) — just a counter.
#
# String representation: Python strings; zip for parallel iteration.
#
# Where this appears in bioinformatics:
# Hamming distance measures genetic distance between sequences of equal
# length — e.g., comparing aligned DNA reads to a reference, counting
# point mutations between homologous genes, or quantifying barcode
# differences in multiplexed sequencing.

import sys
import os


def hamming_distance(p, q):
    """Return the Hamming distance between equal-length strings p and q.

    Args:
        p: First string.
        q: Second string (same length as p).

    Returns:
        Number of positions where p and q differ.
    """
    return sum(1 for a, b in zip(p, q) if a != b)


def main():
    input_file = 'ba1g_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().strip().splitlines()
    else:
        lines = sys.stdin.read().strip().splitlines()

    p = lines[0].strip()
    q = lines[1].strip()

    print(hamming_distance(p, q))


if __name__ == '__main__':
    main()
