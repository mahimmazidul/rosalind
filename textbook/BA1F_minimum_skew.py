# BA1F (Minimum Skew)
# https://rosalind.info/problems/ba1f/
#
# What the algorithm is:
# Find all positions i (0 to |Genome|) where the skew of prefix i is
# minimized. Skew = (#G − #C) in the prefix. Positions of minimum skew
# often indicate the origin of replication in bacterial genomes.
#
# How it works conceptually:
# Initialize skew = 0 at position 0. Scan Genome left to right:
#   - 'G' increases skew by 1
#   - 'C' decreases skew by 1
#   - 'A', 'T' leave skew unchanged
# Record skew at each position. Track the minimum value and all
# positions where it occurs.
#
# Time and space complexity:
# Time: O(n) — single pass over the genome.
# Space: O(n) — to store positions of minimum skew.
#
# String representation: Python string; single character scan.
#
# Where this appears in bioinformatics:
# Minimum skew positions are a classic predictor of the origin of
# replication (oriC) in bacterial genomes — GC skew changes sign near
# the replication origin and termination point.

import sys
import os


def minimum_skew(genome):
    """Return all positions minimizing the G-C skew of each prefix.

    Args:
        genome: A DNA string.

    Returns:
        A list of integer positions (0 to len(genome)) where skew is minimum.
    """
    skew = 0
    min_skew = 0
    positions = [0]  # position 0 has skew 0

    for i, base in enumerate(genome, 1):  # i = 1..n (1-indexed position after reading base i)
        if base == 'G':
            skew += 1
        elif base == 'C':
            skew -= 1

        if skew < min_skew:
            min_skew = skew
            positions = [i]
        elif skew == min_skew:
            positions.append(i)

    return positions


def main():
    input_file = 'ba1f_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().strip().splitlines()
    else:
        lines = sys.stdin.read().strip().splitlines()

    genome = lines[0].strip()

    positions = minimum_skew(genome)
    print(' '.join(str(p) for p in positions))


if __name__ == '__main__':
    main()
