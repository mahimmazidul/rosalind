# PTRA (Protein Translation)
# https://rosalind.info/problems/ptra/
#
# What the genetic code is:
# The genetic code is the set of rules by which information encoded in
# DNA or RNA is translated into proteins. It maps each three-nucleotide
# codon to one of 20 amino acids (or a stop signal). With 4 bases and
# codons of length 3, there are 64 possible codons but only 20 amino
# acids, so the code is redundant — multiple codons can encode the same
# amino acid.
#
# Genetic code variants:
# While the standard genetic code is used by most organisms, some
# organisms and organelles use alternative codes. For example, vertebrate
# mitochondria treat AGA and AGG as stop codons instead of arginine.
# Each variant is assigned a number (1 = standard, 2 = vertebrate
# mitochondrial, etc.). Knowing which code was used is essential for
# correct translation.
#
# What this solution does:
# Given a DNA string and the protein it translates into, this solution
# tests each genetic code variant (tables 1 through 30) using Biopython's
# translate() method, ignoring stop codons, and returns the variant number
# that produces the target protein.


from Bio.Seq import translate
import sys
import os


def find_translation_table(dna, protein):
    """Find the genetic code table that translates DNA to the given protein.

    Tests tables 1 through 30. For each table, translates the DNA string
    starting from the first base, ignoring stop codons, and compares the
    result to the target protein.

    Args:
        dna: A DNA string.
        protein: The expected protein string.

    Returns:
        The table number (integer), or None if no table matches.
    """
    for table in range(1, 31):
        try:
            result = translate(dna, table=table, to_stop=False)
            # Ignore stop codons by removing stop symbols
            cleaned = result.replace('*', '')
            if cleaned == protein:
                return table
        except Exception:
            continue
    return None


def main():
    input_file = 'ptra_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().splitlines()
    else:
        lines = sys.stdin.read().splitlines()

    # Filter out empty lines
    lines = [line.strip() for line in lines if line.strip()]
    dna = lines[0]
    protein = lines[1]

    table = find_translation_table(dna, protein)
    if table is not None:
        print(table)
    else:
        print("No matching table found")


if __name__ == '__main__':
    main()
