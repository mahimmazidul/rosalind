# BA1A (Pattern Count)
# https://rosalind.info/problems/ba1a/
#
# What the algorithm is:
# PatternCount — count how many times a k-mer Pattern appears as a
# substring of Text, including overlapping occurrences.
#
# How it works conceptually:
# Slide a window of length |Pattern| across Text, one position at a
# time, and compare each k-mer substring to Pattern. Count matches.
#
# Time and space complexity:
# Time: O(|Text| * |Pattern|) — naive string comparison at each window.
# Space: O(1) — just a counter.
#
# String representation: Plain Python strings; slicing extracts k-mers.
#
# Where this appears in bioinformatics:
# Pattern counting is the foundation of motif finding — e.g., counting
# occurrences of a transcription factor binding site in a promoter
# region, or finding how often a specific k-mer appears in a genome.

import sys
import os


def pattern_count(text, pattern):
    """Count occurrences of pattern in text, including overlapping.

    Args:
        text: A DNA string.
        pattern: A k-mer pattern to search for.

    Returns:
        The number of times pattern appears in text.
    """
    count = 0
    k = len(pattern)
    for i in range(len(text) - k + 1):
        if text[i:i + k] == pattern:
            count += 1
    return count


def main():
    input_file = 'ba1a_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().strip().splitlines()
    else:
        lines = sys.stdin.read().strip().splitlines()

    text = lines[0].strip()
    pattern = lines[1].strip()

    print(pattern_count(text, pattern))


if __name__ == '__main__':
    main()
