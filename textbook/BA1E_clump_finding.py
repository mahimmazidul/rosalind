# BA1E (Clump Finding)
# https://rosalind.info/problems/ba1e/
#
# What the algorithm is:
# Find all distinct k-mers that form an (L, t)-clump: a k-mer appears
# at least t times within some interval of length L in Genome.
#
# How it works conceptually:
# 1. Collect all starting positions of every k-mer in Genome.
# 2. For each k-mer, sort its positions and slide a window of t
#    consecutive occurrences. If the span from the first to the t-th
#    position fits within L (specifically, pos[j+t-1] - pos[j] <= L - k),
#    the k-mer forms a clump.
# 3. Output all distinct k-mers that form a clump.
#
# Time and space complexity:
# Time: O(|Genome| * k) to collect positions, plus O(total positions)
# for the sliding window checks.
# Space: O(|Genome| * k) — dictionary mapping k-mers to position lists.
#
# String representation: Python strings as dict keys; position lists as values.
#
# Where this appears in bioinformatics:
# Clump finding detects locally overrepresented k-mers — e.g., finding
# putative replication origins (oriC) where specific k-mers cluster
# near the start of replication, or identifying regions of biased
# k-mer composition in genomes.

import sys
import os
from collections import defaultdict


def clump_finding(genome, k, L, t):
    """Find all distinct k-mers forming (L, t)-clumps in genome.

    Args:
        genome: The DNA string to search.
        k: Length of k-mers.
        L: Window length for clump detection.
        t: Minimum occurrences within a window.

    Returns:
        A set of k-mer strings forming clumps.
    """
    positions = defaultdict(list)
    for i in range(len(genome) - k + 1):
        kmer = genome[i:i + k]
        positions[kmer].append(i)

    clumps = set()
    max_span = L - k  # max distance between first and t-th occurrence

    for kmer, pos_list in positions.items():
        if len(pos_list) < t:
            continue
        # pos_list is already sorted (we appended in increasing i)
        for j in range(len(pos_list) - t + 1):
            if pos_list[j + t - 1] - pos_list[j] <= max_span:
                clumps.add(kmer)
                break

    return clumps


def main():
    input_file = 'ba1e_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().strip().splitlines()
    else:
        lines = sys.stdin.read().strip().splitlines()

    genome = lines[0].strip()
    k, L, t = map(int, lines[1].split())

    clumps = clump_finding(genome, k, L, t)
    print(' '.join(sorted(clumps)))


if __name__ == '__main__':
    main()
