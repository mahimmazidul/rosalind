# BA1I (Frequent Words with Mismatches)
# https://rosalind.info/problems/ba1i/
#
# What the algorithm is:
# Find all k-mers that maximize Count_d(Text, Pattern) — the number of
# occurrences of Pattern in Text with at most d mismatches. Pattern
# need not appear exactly in Text.
#
# How it works conceptually:
# For each k-mer substring of Text, generate its entire d-neighborhood
# (all k-mers within Hamming distance d). Tally each neighbor across
# all windows. The k-mers with the highest tally are the answer.
#
# Why neighborhoods instead of all 4^k k-mers: for k ≤ 12, d ≤ 3, the
# neighborhood per window is small (~6K for k=12,d=3), far fewer than
# 4^12 ≈ 16M.
#
# Time and space complexity:
# Time: O(|Text| * neighborhood_size * k) — generate neighborhoods for
# each window and tally.
# Space: O(number of distinct neighbors) — dictionary of counts.
#
# String representation: Python strings; recursive neighborhood generation.
#
# Where this appears in bioinformatics:
# Finding frequent k-mers with mismatches discovers regulatory motifs
# that occur with slight variations across promoters — robust to SNPs
# and sequencing errors.

import sys
import os
from collections import defaultdict
from itertools import product


def hamming_distance(p, q):
    """Return the Hamming distance between equal-length strings p and q."""
    return sum(1 for a, b in zip(p, q) if a != b)


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

    # Case 1: keep first character, vary suffix with up to d mismatches
    for suffix_neighbor in neighbors(suffix, d):
        yield first + suffix_neighbor

    # Case 2: change first character (uses 1 mismatch), vary suffix with up to d-1
    if d > 0:
        for base in 'ACGT':
            if base != first:
                for suffix_neighbor in neighbors(suffix, d - 1):
                    yield base + suffix_neighbor


def frequent_words_with_mismatches(text, k, d):
    """Return all most frequent k-mers with up to d mismatches in text.

    Args:
        text: A DNA string.
        k: Length of k-mers.
        d: Maximum mismatches per occurrence.

    Returns:
        A list of k-mer strings with the highest Count_d value.
    """
    counts = defaultdict(int)

    for i in range(len(text) - k + 1):
        window = text[i:i + k]
        for neighbor in neighbors(window, d):
            counts[neighbor] += 1

    max_count = max(counts.values())
    return [kmer for kmer, cnt in counts.items() if cnt == max_count]


def main():
    input_file = 'ba1i_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().strip().splitlines()
    else:
        lines = sys.stdin.read().strip().splitlines()

    text = lines[0].strip()
    k, d = map(int, lines[1].split())

    result = frequent_words_with_mismatches(text, k, d)
    print(' '.join(sorted(result)))


if __name__ == '__main__':
    main()
