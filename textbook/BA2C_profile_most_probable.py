# BA2C (Profile-most Probable k-mer)
# https://rosalind.info/problems/ba2c/
#
# What the algorithm is:
# Find the k-mer in Text with the highest probability under a given
# 4×k profile matrix (rows: A, C, G, T; columns: positions).
#
# How it works conceptually:
# For each k-mer window of Text, compute its probability as the product
# of Profile[symbol][position] across all positions. Return the window
# with the highest probability (first one in case of ties).
#
# Time and space complexity:
# Time: O(|Text| * k) — evaluate each window.
# Space: O(k) — the profile matrix and best k-mer.
#
# String representation: Symbol-to-row index mapping (A=0,C=1,G=2,T=3).
#
# Where this appears in bioinformatics:
# Profile-most probable k-mer is the core scoring step of Gibbs
# sampling and EM-based motif discovery — finding the best motif
# instance in each sequence given a current position weight matrix.

import sys
import os


def profile_most_probable(text, k, profile):
    """Return the Profile-most probable k-mer in text.

    Args:
        text: A DNA string.
        k: Length of k-mer.
        profile: A 4×k list of lists; rows are A,C,G,T probabilities.

    Returns:
        The k-mer string with the highest profile probability.
    """
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


def main():
    input_file = 'ba2c_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().strip().splitlines()
    else:
        lines = sys.stdin.read().strip().splitlines()

    text = lines[0].strip()
    k = int(lines[1].strip())
    profile = []
    for i in range(4):
        profile.append(list(map(float, lines[2 + i].split())))

    result = profile_most_probable(text, k, profile)
    print(result)


if __name__ == '__main__':
    main()
