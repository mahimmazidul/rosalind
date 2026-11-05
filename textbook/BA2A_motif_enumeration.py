# BA2A (Motif Enumeration)
# https://rosalind.info/problems/ba2a/
#
# What the algorithm is:
# Find all (k, d)-motifs: k-mers that appear in every string of Dna
# with at most d mismatches. For each k-mer window in Dna, generate
# its d-neighborhood and test each candidate against all strings.
#
# How it works conceptually:
# 1. Collect all k-mer windows from all strings in Dna.
# 2. For each window, generate its d-neighborhood (all k-mers within d
#    mismatches).
# 3. For each candidate neighbor, check if it appears in every string
#    of Dna with at most d mismatches (i.e., some k-mer window in each
#    string is within d of the candidate).
# 4. Collect all passing candidates, remove duplicates.
#
# Time and space complexity:
# Time: O(|Dna| * n * neighborhood_size * |Dna| * n * k) — worst case,
# but practical for small k, d.
# Space: O(neighborhood_size) — set of candidate motifs.
#
# String representation: Hamming distance for approximate matching.
#
# Where this appears in bioinformatics:
# Implanted motif discovery — finding a short regulatory motif (e.g. a
# transcription factor binding site) that occurs with slight variations
# across multiple promoter sequences.

import sys
import os
from itertools import product


def hamming_distance(p, q):
    """Return the Hamming distance between equal-length strings p and q."""
    return sum(1 for a, b in zip(p, q) if a != b)


def appears_with_mismatches(pattern, text, d):
    """Check if pattern appears in text with at most d mismatches."""
    k = len(pattern)
    for i in range(len(text) - k + 1):
        if hamming_distance(pattern, text[i:i + k]) <= d:
            return True
    return False


def neighbors(pattern, d):
    """Generate all k-mers within Hamming distance d of pattern."""
    if d == 0:
        yield pattern
        return
    if len(pattern) == 0:
        yield ''
        return

    first = pattern[0]
    suffix = pattern[1:]

    for suffix_neighbor in neighbors(suffix, d):
        yield first + suffix_neighbor

    if d > 0:
        for base in 'ACGT':
            if base != first:
                for suffix_neighbor in neighbors(suffix, d - 1):
                    yield base + suffix_neighbor


def motif_enumeration(dna, k, d):
    """Find all (k, d)-motifs in a collection of DNA strings.

    Args:
        dna: A list of DNA strings.
        k: Length of the motif.
        d: Maximum mismatches allowed.

    Returns:
        A set of k-mer strings that are (k, d)-motifs.
    """
    motifs = set()

    # Collect all k-mer windows from all strings
    windows = []
    for text in dna:
        for i in range(len(text) - k + 1):
            windows.append(text[i:i + k])

    # For each window, generate its d-neighborhood and test candidates
    seen = set()
    for window in windows:
        for candidate in neighbors(window, d):
            if candidate in seen:
                continue
            seen.add(candidate)
            if all(appears_with_mismatches(candidate, text, d) for text in dna):
                motifs.add(candidate)

    return motifs


def main():
    input_file = 'ba2a_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().strip().splitlines()
    else:
        lines = sys.stdin.read().strip().splitlines()

    k, d = map(int, lines[0].split())
    dna = [line.strip() for line in lines[1:] if line.strip()]

    motifs = motif_enumeration(dna, k, d)
    print(' '.join(sorted(motifs)))


if __name__ == '__main__':
    main()
