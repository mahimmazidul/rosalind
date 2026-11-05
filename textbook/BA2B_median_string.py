# BA2B (Median String)
# https://rosalind.info/problems/ba2b/
#
# What the algorithm is:
# Find a k-mer (median string) that minimizes the total distance to a
# collection of DNA strings. Distance from a pattern to a string is the
# minimum Hamming distance between the pattern and any k-mer in that
# string; the total distance sums this over all strings.
#
# How it works conceptually:
# 1. Generate all 4^k possible k-mers (A=0,C=1,G=2,T=3 in base-4).
# 2. For each k-mer, compute its distance to each DNA string (min Hamming
#    distance to any k-mer window in that string), then sum.
# 3. Return the k-mer with the smallest total distance.
#
# Time and space complexity:
# Time: O(4^k * |Dna| * n * k) — enumerate all k-mers, scan each string.
# Space: O(k) — best pattern and its score.
#
# String representation: Base-4 enumeration via NumberToPattern.
#
# Where this appears in bioinformatics:
# The median string is a principled motif-finding approach — finding the
# k-mer that best represents a set of sequences, robust to variations
# at each position.

import sys
import os
from itertools import product


def hamming_distance(p, q):
    """Return the Hamming distance between equal-length strings p and q."""
    return sum(1 for a, b in zip(p, q) if a != b)


def distance_pattern_to_string(pattern, text):
    """Minimum Hamming distance between pattern and any k-mer in text."""
    k = len(pattern)
    min_dist = k  # max possible distance
    for i in range(len(text) - k + 1):
        dist = hamming_distance(pattern, text[i:i + k])
        if dist < min_dist:
            min_dist = dist
    return min_dist


def distance_pattern_to_dna(pattern, dna):
    """Sum of distances from pattern to each string in dna."""
    return sum(distance_pattern_to_string(pattern, text) for text in dna)


def number_to_pattern(index, k):
    """Convert an integer index to a DNA string of length k."""
    num_to_symbol = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    digits = []
    for _ in range(k):
        digits.append(num_to_symbol[index % 4])
        index //= 4
    return ''.join(reversed(digits))


def median_string(dna, k):
    """Find a k-mer minimizing total distance to dna.

    Args:
        dna: A list of DNA strings.
        k: Length of the motif.

    Returns:
        A k-mer string that minimizes d(Pattern, Dna).
    """
    best_pattern = None
    best_distance = float('inf')

    for i in range(4 ** k):
        pattern = number_to_pattern(i, k)
        dist = distance_pattern_to_dna(pattern, dna)
        if dist < best_distance:
            best_distance = dist
            best_pattern = pattern

    return best_pattern


def main():
    input_file = 'ba2b_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().strip().splitlines()
    else:
        lines = sys.stdin.read().strip().splitlines()

    k = int(lines[0].strip())
    dna = [line.strip() for line in lines[1:] if line.strip()]

    result = median_string(dna, k)
    print(result)


if __name__ == '__main__':
    main()
