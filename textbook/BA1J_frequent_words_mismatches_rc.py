# BA1J (Frequent Words with Mismatches and Reverse Complements)
# https://rosalind.info/problems/ba1j/
#
# What the algorithm is:
# Find all k-mers Pattern maximizing Count_d(Text, Pattern) +
# Count_d(Text, reverse_complement(Pattern)). A k-mer and its reverse
# complement are treated as equivalent signals since DNA is double-stranded.
#
# How it works conceptually:
# For each k-mer window W of Text, generate its d-neighborhood N(W,d)
# (all k-mers within Hamming distance d of W). For each P in N(W,d):
#   - score[P] += 1            (P matches W in Text)
#   - score[rc(P)] += 1        (rc(P) matches W in Text, i.e. P matches rc(W))
# The second term effectively counts matches against rc(Text), since
# HammingDistance(rc(P), W) = HammingDistance(P, rc(W)).
# Finally, output all k-mers with maximum score.
#
# Time and space complexity:
# Time: O(|Text| * neighborhood_size * k)
# Space: O(distinct neighbors) — dictionary of scores.
#
# String representation: Python strings; reverse complement via lookup dict.
#
# Where this appears in bioinformatics:
# Motif finding on double-stranded DNA — a transcription factor binding
# site on one strand implies its reverse complement on the other, so
# both strands are searched simultaneously.

import sys
import os
from collections import defaultdict


def hamming_distance(p, q):
    """Return the Hamming distance between equal-length strings p and q."""
    return sum(1 for a, b in zip(p, q) if a != b)


def reverse_complement(pattern):
    """Return the reverse complement of a DNA string."""
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(comp[b] for b in reversed(pattern))


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


def frequent_words_mismatches_rc(text, k, d):
    """Return all k-mers maximizing Count_d(text, P) + Count_d(text, rc(P))."""
    score = defaultdict(int)

    for i in range(len(text) - k + 1):
        window = text[i:i + k]
        for p in neighbors(window, d):
            score[p] += 1
            score[reverse_complement(p)] += 1

    max_score = max(score.values())
    return [kmer for kmer, s in score.items() if s == max_score]


def main():
    input_file = 'ba1j_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().strip().splitlines()
    else:
        lines = sys.stdin.read().strip().splitlines()

    text = lines[0].strip()
    k, d = map(int, lines[1].split())

    result = frequent_words_mismatches_rc(text, k, d)
    print(' '.join(sorted(result)))


if __name__ == '__main__':
    main()
