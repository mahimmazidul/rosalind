# BA1B (Most Frequent Words / k-mers)
# https://rosalind.info/problems/ba1b/
#
# What the algorithm is:
# Find all most frequent k-mers in a DNA string — every k-mer whose
# occurrence count equals the maximum count among all k-mers.
#
# How it works conceptually:
# Slide a window of length k across Text, tally each k-mer in a
# dictionary. Find the maximum count, then output all k-mers that
# reach it.
#
# Time and space complexity:
# Time: O(|Text| * k) — each of O(|Text|) k-mers takes O(k) to hash.
# Space: O(|Text| * k) — dictionary of distinct k-mers.
#
# String representation: Python strings as dictionary keys.
#
# Where this appears in bioinformatics:
# Finding frequent k-mers is the first step in motif discovery —
# discovering overrepresented sequence patterns (e.g. transcription
# factor binding sites, replication origins) in genomic DNA.

import sys
import os


def frequent_words(text, k):
    """Return all most frequent k-mers in text (in any order).

    Args:
        text: A DNA string.
        k: Length of k-mers to search for.

    Returns:
        A list of k-mer strings with the highest occurrence count.
    """
    counts = {}
    for i in range(len(text) - k + 1):
        kmer = text[i:i + k]
        counts[kmer] = counts.get(kmer, 0) + 1

    max_count = max(counts.values())
    return [kmer for kmer, cnt in counts.items() if cnt == max_count]


def main():
    input_file = 'ba1b_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().strip().splitlines()
    else:
        lines = sys.stdin.read().strip().splitlines()

    text = lines[0].strip()
    k = int(lines[1].strip())

    result = frequent_words(text, k)
    print(' '.join(result))


if __name__ == '__main__':
    main()
