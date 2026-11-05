# BA1K (Frequency Array)
# https://rosalind.info/problems/ba1k/
#
# What the algorithm is:
# Generate the frequency array of k-mers in a DNA string — an array of
# length 4^k where index i holds the count of the i-th k-mer in
# lexicographic order (A < C < G < T).
#
# How it works conceptually:
# Map each k-mer to a base-4 integer (A=0, C=1, G=2, T=3), treating
# the string as a base-4 number. Use this integer as an index into a
# frequency array. Slide a k-mer window across Text and increment the
# corresponding index.
#
# Time and space complexity:
# Time: O(|Text| * k) — each window takes O(k) to convert to a number.
# Space: O(4^k) — the frequency array.
#
# String representation: Base-4 encoding of k-mers as array indices.
#
# Where this appears in bioinformatics:
# Frequency arrays are the compact representation behind k-mer spectra
# — used in genome assembly (k-mer frequency histograms), metagenomic
# binning, and sequence composition analysis.

import sys
import os


def pattern_to_number(pattern):
    """Convert a DNA k-mer to its base-4 index (lexicographic order).

    Args:
        pattern: A DNA string over {A, C, G, T}.

    Returns:
        An integer in [0, 4^k - 1].
    """
    symbol_to_num = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    number = 0
    for base in pattern:
        number = number * 4 + symbol_to_num[base]
    return number


def frequency_array(text, k):
    """Return the frequency array of k-mers in text.

    Args:
        text: A DNA string.
        k: Length of k-mers.

    Returns:
        A list of length 4^k with k-mer counts in lexicographic order.
    """
    size = 4 ** k
    freq = [0] * size
    for i in range(len(text) - k + 1):
        kmer = text[i:i + k]
        idx = pattern_to_number(kmer)
        freq[idx] += 1
    return freq


def main():
    input_file = 'ba1k_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().strip().splitlines()
    else:
        lines = sys.stdin.read().strip().splitlines()

    text = lines[0].strip()
    k = int(lines[1].strip())

    freq = frequency_array(text, k)
    print(' '.join(str(x) for x in freq))


if __name__ == '__main__':
    main()
