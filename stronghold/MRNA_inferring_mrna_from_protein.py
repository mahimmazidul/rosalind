# MRNA (Inferring mRNA from Protein)
# https://rosalind.info/problems/mrna/
#
# What this problem is about:
# Translation is the process where a ribosome reads mRNA in triplets called
# codons and builds a protein from amino acids. Going forward (RNA → protein)
# is deterministic: each codon maps to exactly one amino acid (or stop).
# Going backward (protein → RNA) is not: most amino acids are encoded by
# multiple codons. Leucine has 6 codons, serine has 6, arginine has 6,
# while methionine and tryptophan each have just 1. This redundancy is
# called the degeneracy of the genetic code.
#
# Why this matters:
# When researchers find a new protein and want to find the gene that
# encodes it, they face a huge search space. A protein of just 100 amino
# acids could be encoded by an astronomical number of different RNA strings
# because of codon degeneracy. This is why we can't simply "reverse"
# translation — we need to count the possibilities and work with modular
# arithmetic to keep the numbers manageable.
#
# What this solution does:
# For each amino acid in the protein, count how many codons encode it.
# Multiply all those counts together, then multiply by 3 (the number of
# stop codons: UAA, UAG, UGA — every translated RNA must end with one).
# Take the result modulo 1,000,000 to keep the number manageable.


# Number of codons encoding each amino acid, from the standard genetic code.
CODON_COUNTS = {
    'A': 4,  # GCU, GCC, GCA, GCG
    'R': 6,  # CGU, CGC, CGA, CGG, AGA, AGG
    'N': 2,  # AAU, AAC
    'D': 2,  # GAU, GAC
    'C': 2,  # UGU, UGC
    'Q': 2,  # CAA, CAG
    'E': 2,  # GAA, GAG
    'G': 4,  # GGU, GGC, GGA, GGG
    'H': 2,  # CAU, CAC
    'I': 3,  # AUU, AUC, AUA
    'L': 6,  # UUA, UUG, CUU, CUC, CUA, CUG
    'K': 2,  # AAA, AAG
    'M': 1,  # AUG
    'F': 2,  # UUU, UUC
    'P': 4,  # CCU, CCC, CCA, CCG
    'S': 6,  # UCU, UCC, UCA, UCG, AGU, AGC
    'T': 4,  # ACU, ACC, ACA, ACG
    'W': 1,  # UGG
    'Y': 2,  # UAU, UAC
    'V': 4,  # GUU, GUC, GUA, GUG
}

STOP_CODONS = 3  # UAA, UAG, UGA
MODULO = 1_000_000


def count_rna_strings(protein):
    """Return the number of RNA strings that translate to the given protein,
    modulo 1,000,000.

    Each amino acid contributes its codon count. The RNA must end with a
    stop codon, contributing a factor of 3.

    Args:
        protein: A string of amino acid single-letter codes.

    Returns:
        The count modulo 1,000,000.
    """
    result = STOP_CODONS
    for aa in protein:
        result = (result * CODON_COUNTS[aa]) % MODULO
    return result


def main():
    import sys
    import os

    input_file = 'mrna_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            protein = f.read().strip()
    else:
        protein = sys.stdin.read().strip()

    print(count_rna_strings(protein))


if __name__ == '__main__':
    main()
