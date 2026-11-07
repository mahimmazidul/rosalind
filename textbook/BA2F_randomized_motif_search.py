# BA2F (Randomized Motif Search)
# https://rosalind.info/problems/ba2f/
#
# What the algorithm is:
# Randomized Motif Search — a Monte Carlo algorithm for motif finding.
# Repeatedly starts from random k-mers, iteratively refines via profile
# and scoring, and keeps the best result across 1000 runs.
#
# How it works conceptually:
# Each run:
# 1. Randomly pick one k-mer from each DNA string as initial motifs.
# 2. Loop:
#    a. Build a profile matrix from current motifs (with pseudocounts).
#    b. Replace each motif with the profile-most probable k-mer in its
#       string.
#    c. If the new motif set scores better, keep it; otherwise stop.
# 3. Track the best motif set across all 1000 runs.
#
# Score: sum over columns of (t - count of most frequent symbol).
# Lower = more conserved.
#
# Time and space complexity:
# Time: O(1000 * t * n * k * 4) per run in the worst case.
# Space: O(t * k) — motif matrix.
#
# String representation: Profile matrix (4×k) with pseudocounts.
#
# Where this appears in bioinformatics:
# Randomized motif search is used for de novo motif discovery —
# finding shared regulatory sequences across a set of unaligned
# promoter regions without knowing the motif in advance.

import sys
import os
import random


def profile_most_probable(text, k, profile):
    """Return the first k-mer in text with the highest profile probability."""
    symbol_to_row = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    best_kmer = None
    best_prob = -1.0
    for i in range(len(text) - k + 1):
        kmer = text[i:i + k]
        prob = 1.0
        for j, base in enumerate(kmer):
            prob *= profile[symbol_to_row[base]][j]
        if prob > best_prob:
            best_prob = prob
            best_kmer = kmer
    return best_kmer


def build_profile(motifs):
    """Build a 4×k profile matrix with pseudocounts from a list of motifs."""
    t = len(motifs)
    k = len(motifs[0])
    symbol_to_row = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    counts = [[1] * k for _ in range(4)]

    for motif in motifs:
        for j, base in enumerate(motif):
            counts[symbol_to_row[base]][j] += 1

    profile = [[counts[r][j] / (t + 4) for j in range(k)] for r in range(4)]
    return profile


def score_motifs(motifs):
    """Score a set of motifs: sum over columns of (t - max_count)."""
    t = len(motifs)
    k = len(motifs[0])
    score = 0
    for j in range(k):
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for motif in motifs:
            counts[motif[j]] += 1
        max_count = max(counts.values())
        score += t - max_count
    return score


def randomized_motif_search(dna, k, t):
    """Single run of Randomized Motif Search.

    Args:
        dna: A list of t DNA strings.
        k: Length of the motif.
        t: Number of strings.

    Returns:
        A list of t k-mer strings (BestMotifs for this run).
    """
    # Randomly select k-mers from each string
    motifs = []
    for text in dna:
        i = random.randint(0, len(text) - k)
        motifs.append(text[i:i + k])

    best_motifs = motifs[:]
    best_score = score_motifs(motifs)

    while True:
        profile = build_profile(motifs)
        motifs = [profile_most_probable(text, k, profile) for text in dna]
        current_score = score_motifs(motifs)
        if current_score < best_score:
            best_score = current_score
            best_motifs = motifs[:]
        else:
            return best_motifs


def main():
    input_file = 'ba2f_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().strip().splitlines()
    else:
        lines = sys.stdin.read().strip().splitlines()

    k, t = map(int, lines[0].split())
    dna = [line.strip() for line in lines[1:] if line.strip()]

    best_overall = None
    best_overall_score = float('inf')

    for _ in range(1000):
        motifs = randomized_motif_search(dna, k, t)
        score = score_motifs(motifs)
        if score < best_overall_score:
            best_overall_score = score
            best_overall = motifs[:]

    print('\n'.join(best_overall))


if __name__ == '__main__':
    main()
