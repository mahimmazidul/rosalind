# BA1D (Pattern Matching)
# https://rosalind.info/problems/ba1d/
#
# What the algorithm is:
# Find all starting positions (0-based) where Pattern occurs as a
# substring of Genome, including overlapping occurrences.
#
# How it works conceptually:
# Slide a window of length |Pattern| across Genome. At each position,
# check if the substring equals Pattern. Record every match position.
#
# Time and space complexity:
# Time: O(|Genome| * |Pattern|) — naive comparison at each window.
# Space: O(|Genome|) — list of match positions.
#
# String representation: Python strings; slicing extracts k-mers.
#
# Where this appears in bioinformatics:
# Pattern matching locates a motif (e.g. a restriction site or binding
# site) across a genome — reporting every position where it occurs,
# including overlapping hits on both strands.

import sys
import os


def pattern_matching(pattern, genome):
    """Return all 0-based start positions where pattern occurs in genome.

    Args:
        pattern: The pattern string to search for.
        genome: The text to search in.

    Returns:
        A list of integer starting positions (0-based).
    """
    positions = []
    k = len(pattern)
    for i in range(len(genome) - k + 1):
        if genome[i:i + k] == pattern:
            positions.append(i)
    return positions


def main():
    input_file = 'ba1d_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().strip().splitlines()
    else:
        lines = sys.stdin.read().strip().splitlines()

    pattern = lines[0].strip()
    genome = lines[1].strip()

    positions = pattern_matching(pattern, genome)
    print(' '.join(str(p) for p in positions))


if __name__ == '__main__':
    main()
