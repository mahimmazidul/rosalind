# INI (Introduction to the Bioinformatics Armory)
# https://rosalind.info/problems/ini/
#
# What the tool is:
# The Sequence Manipulation Suite (SMS 2) is a collection of web-based
# programs for generating, formatting, and analyzing short DNA and protein
# sequences. One of its simplest tools, "DNA stats," counts how many times
# each nucleotide (A, C, G, T) appears in a DNA strand. In the Armory track,
# we use BioPython as the Python equivalent of SMS 2 — it's a free, open
# source library that provides biologically aware data structures like Seq,
# which extends Python strings with methods such as count(), translate(),
# and reverse_complement().
#
# Why this tool exists in real bioinformatics workflows:
# Biologists constantly need to count nucleotides — to check sequencing
# quality, compute GC content, or spot biases in a sequencing run. Rather
# than writing a custom script every time, tools like SMS 2 and BioPython
# provide reliable, tested functions that work out of the box. BioPython
# in particular lets Python programmers handle biological data with
# domain-appropriate data structures instead of treating DNA as plain text.
#
# What this solution does:
# It reads a DNA string and counts the occurrences of A, C, G, and T using
# BioPython's Seq object, then prints the four counts separated by spaces.


from Bio.Seq import Seq


def count_nucleotides(dna):
    """Count A, C, G, T nucleotides in a DNA string using BioPython's Seq.

    Args:
        dna: A string of 'A', 'C', 'G', 'T' characters.

    Returns:
        A tuple of (count_A, count_C, count_G, count_T).
    """
    seq = Seq(dna)
    return seq.count('A'), seq.count('C'), seq.count('G'), seq.count('T')


def main():
    import sys
    import os

    input_file = 'ini_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            dna = f.read().strip()
    else:
        dna = sys.stdin.read().strip()

    a, c, g, t = count_nucleotides(dna)
    print(f'{a} {c} {g} {t}')


if __name__ == '__main__':
    main()
