# BA1H (Approximate Pattern Matching)
# https://rosalind.info/problems/ba1h/
#
# What the algorithm is:
# Find all starting positions where Pattern appears as a substring of
# Text with at most d mismatches (Hamming distance ≤ d).
#
# How it works conceptually:
# Slide a window of length |Pattern| across Text. At each position,
# compute the Hamming distance between Pattern and the window. Record
# positions where the distance is ≤ d.
#
# Time and space complexity:
# Time: O(|Text| * |Pattern|) — Hamming distance at each window.
# Space: O(|Text|) — list of match positions.
#
# String representation: Python strings; slicing extracts windows.
#
# Where this appears in bioinformatics:
# Approximate pattern matching finds binding sites or motifs that
# occur with minor variations (SNPs, sequencing errors) — e.g., locating
# a transcription factor binding site that may differ by a few bases.

import sys
import os


def hamming_distance(p, q):
    """Return the Hamming distance between equal-length strings p and q."""
    return sum(1 for a, b in zip(p, q) if a != b)


def approximate_pattern_matching(pattern, text, d):
    """Return all start positions where pattern matches text with ≤ d mismatches.

    Args:
        pattern: The pattern string.
        text: The text to search in.
        d: Maximum allowed mismatches.

    Returns:
        A list of 0-based starting positions.
    """
    positions = []
    k = len(pattern)
    for i in range(len(text) - k + 1):
        if hamming_distance(pattern, text[i:i + k]) <= d:
            positions.append(i)
    return positions


def main():
    input_file = 'ba1h_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().strip().splitlines()
    else:
        lines = sys.stdin.read().strip().splitlines()

    pattern = lines[0].strip()
    text = lines[1].strip()
    d = int(lines[2].strip())

    positions = approximate_pattern_matching(pattern, text, d)
    print(' '.join(str(p) for p in positions))


if __name__ == '__main__':
    main()
