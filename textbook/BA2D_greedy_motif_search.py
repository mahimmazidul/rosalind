# BA2D (Greedy Motif Search)
# https://rosalind.info/problems/ba2d/
#
# What the algorithm is:
# Greedy Motif Search finds a set of motifs (one k-mer per DNA string)
# that minimizes the motif score. It starts from every possible k-mer
# in the first string, builds a profile from the motifs found so far
# (with pseudocounts), picks the profile-most probable k-mer in each
# subsequent string, and keeps the best-scoring set.
#
# How it works conceptually:
# 1. Initialize BestMotifs with the first k-mer from each string.
# 2. For each k-mer in the first string:
#    a. Start Motifs with that k-mer.
#    b. For each remaining string, build a profile from current motifs
#       (no pseudocounts), select the profile-most probable k-mer.
#    c. Score the resulting motif set; keep it if better than BestMotifs.
# 3. Return BestMotifs.
#
# Score: sum over columns of (t - count of most frequent symbol).
# Lower score = more conserved motifs = better.
#
# Time and space complexity:
# Time: O(n * t * n * k * 4) — for each of O(n) starting k-mers, build
# profiles and score for t strings.
# Space: O(t * k) — motif matrix.
#
# String representation: Profile matrix (4×k) with pseudocounts.
#
# Where this appears in bioinformatics:
# Greedy motif search is a heuristic for finding transcription factor
# binding sites — a position weight matrix is built incrementally from
# candidate motif instances across promoter sequences.

import sys
import os


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
    """Build a 4×k profile matrix from a list of motifs (no pseudocounts).

    Args:
        motifs: A list of equal-length DNA strings.

    Returns:
        A 4×k list of lists: rows are A,C,G,T probabilities.
    """
    t = len(motifs)
    k = len(motifs[0])
    symbol_to_row = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    counts = [[0] * k for _ in range(4)]

    for motif in motifs:
        for j, base in enumerate(motif):
            counts[symbol_to_row[base]][j] += 1

    profile = [[counts[r][j] / t for j in range(k)] for r in range(4)]
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


def greedy_motif_search(dna, k, t):
    """Run Greedy Motif Search.

    Args:
        dna: A list of t DNA strings.
        k: Length of the motif.
        t: Number of strings.

    Returns:
        A list of t k-mer strings (BestMotifs).
    """
    # Initialize BestMotifs with first k-mer from each string
    best_motifs = [text[:k] for text in dna]
    best_score = score_motifs(best_motifs)

    # Try each k-mer in the first string as the starting motif
    first_string = dna[0]
    for i in range(len(first_string) - k + 1):
        motifs = [first_string[i:i + k]]
        for j in range(1, t):
            profile = build_profile(motifs)
            motif = profile_most_probable(dna[j], k, profile)
            motifs.append(motif)
        current_score = score_motifs(motifs)
        if current_score < best_score:
            best_score = current_score
            best_motifs = motifs[:]

    return best_motifs


def main():
    input_file = 'ba2d_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().strip().splitlines()
    else:
        lines = sys.stdin.read().strip().splitlines()

    k, t = map(int, lines[0].split())
    dna = [line.strip() for line in lines[1:] if line.strip()]

    motifs = greedy_motif_search(dna, k, t)
    print('\n'.join(motifs))


if __name__ == '__main__':
    main()
