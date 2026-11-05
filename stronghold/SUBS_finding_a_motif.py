# SUBS (Finding a Motif in DNA)
# https://rosalind.info/problems/subs/
#
# What a motif is:
# A motif is a short, recurring interval of DNA that tends to appear in the
# genomes of different organisms or in multiple places within the same genome.
# When the same stretch of DNA shows up in two species, it's a strong hint
# that the sequence has an important biological function — natural selection
# tends to preserve useful sequences while letting junk DNA drift and mutate.
# Motifs can mark binding sites for proteins, regulatory regions that control
# when genes turn on or off, or functional parts of RNAs.
#
# Why finding motifs matters in real biology:
# Genomes are full of repeated sequences — some functional, some parasitic.
# The most common repeat in humans is the Alu element, about 300 base pairs
# long and appearing around a million times in every human genome. Alu hasn't
# been found to serve a positive purpose; in fact, when a new Alu repeat
# inserts itself into a genome, it can cause genetic disorders. Being able to
# scan a genome for a known motif is how researchers find where a regulatory
# protein might bind, where a virus might have inserted itself, or where a
# disease-associated variant might be hiding.
#
# What this solution does:
# Given two DNA strings s and t, it finds every starting position where t
# occurs as a contiguous substring of s. Positions are reported using 1-based
# numbering, the convention used in biology (the first base of a sequence is
# position 1, not position 0 as in Python). This is the kind of "combing
# through the haystack" that molecular biologists do whenever they search a
# genome for a known sequence motif.


def find_motif_positions(s, t):
    """Return all 1-based starting positions where t occurs in s as a substring.

    Args:
        s: The DNA string to search within.
        t: The motif (substring) to search for.

    Returns:
        A list of integers, each a 1-based starting position of t in s.
    """
    positions = []
    t_len = len(t)
    for i in range(len(s) - t_len + 1):
        if s[i:i + t_len] == t:
            positions.append(i + 1)  # convert 0-based to 1-based
    return positions


def main():
    import sys
    import os

    input_file = 'subs_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().splitlines()
    else:
        lines = sys.stdin.read().splitlines()

    # Filter out empty lines and take the first two non-empty lines
    dna = None
    motif = None
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if dna is None:
            dna = line
        elif motif is None:
            motif = line
            break

    positions = find_motif_positions(dna, motif)
    print(' '.join(str(p) for p in positions))


if __name__ == '__main__':
    main()
