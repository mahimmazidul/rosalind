# PROT (Translating RNA into Protein)
# https://rosalind.info/problems/prot/
#
# What translation is:
# Translation is the process where a ribosome reads a messenger RNA (mRNA)
# molecule and builds a protein from it. The ribosome moves along the mRNA
# in steps of three nucleotides. Each three-letter group is called a codon,
# and each codon specifies one amino acid. Transfer RNA (tRNA) molecules
# bring the correct amino acid to the ribosome, matching their anticodon to
# the mRNA codon. The amino acids are joined together into a growing chain
# — a polypeptide that folds into a protein.
#
# Why the genetic code is a triplet code:
# There are four RNA bases (A, U, C, G) and twenty standard amino acids.
# If each codon were one base long, we could only encode 4 amino acids.
# If two bases long, 4^2 = 16 — still not enough. With three bases per
# codon, we get 4^3 = 64 possible codons, which is more than enough for 20
# amino acids. The extra codons mean the code is redundant: several codons
# can encode the same amino acid. For example, leucine is encoded by six
# different codons. This redundancy is called degeneracy, and it buffers
# against mutations — a change in the third base of a codon often doesn't
# change the amino acid.
#
# Start and stop codons:
# Translation begins at a start codon, AUG, which also codes for Methionine
# (M). It's the signal that tells the ribosome where to start reading.
# Translation ends at a stop codon: UAA, UAG, or UGA. These don't code for
# any amino acid — instead, release factors bind to them and the ribosome
# releases the finished polypeptide. In this problem, we translate from the
# start until we hit a stop codon, just like a real ribosome does.
#
# Why the genetic code is nearly universal:
# Almost all organisms on Earth use the same codon table — from bacteria to
# humans. A human gene inserted into bacteria will usually produce the same
# protein. This near-universality is strong evidence that all life shares a
# common ancestor; the code was established early in evolution and has been
# conserved ever since because changing it would break every protein at once.
#
# Practical importance:
# - Recombinant protein production: inserting a human gene into bacteria to
#   make insulin, growth hormone, or vaccines relies on the universal code.
# - Synthetic gene design: scientists can write DNA sequences knowing
#   exactly what protein they'll produce.
# - Understanding mutations: a nucleotide change might be silent, or it
#   might change an amino acid and alter protein function.
#
# This is Rosalind problem #7. Given an RNA string, translate it into a
# protein string using the standard genetic code, stopping at the first stop
# codon. The codon table is written out inline (not imported) so the file is
# self-contained and educational — you can see the full code at a glance.


# Standard genetic code: maps each RNA codon to its amino acid single-letter
# code, or 'Stop' for termination codons. Written inline so this file is
# self-contained and the full code is visible without external imports.
CODON_TABLE = {
    # Alanine
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    # Arginine
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    # Asparagine
    'AAU': 'N', 'AAC': 'N',
    # Aspartic acid
    'GAU': 'D', 'GAC': 'D',
    # Cysteine
    'UGU': 'C', 'UGC': 'C',
    # Glutamine
    'CAA': 'Q', 'CAG': 'Q',
    # Glutamic acid
    'GAA': 'E', 'GAG': 'E',
    # Glycine
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    # Histidine
    'CAU': 'H', 'CAC': 'H',
    # Isoleucine
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I',
    # Leucine
    'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    # Lysine
    'AAA': 'K', 'AAG': 'K',
    # Methionine (start codon)
    'AUG': 'M',
    # Phenylalanine
    'UUU': 'F', 'UUC': 'F',
    # Proline
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    # Serine
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'AGU': 'S', 'AGC': 'S',
    # Threonine
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    # Tryptophan
    'UGG': 'W',
    # Tyrosine
    'UAU': 'Y', 'UAC': 'Y',
    # Valine
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    # Stop codons
    'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop',
}


def translate(rna):
    """Translate an RNA string into a protein string.

    Reads the RNA in codons (groups of 3 bases) and looks up each codon in
    the genetic code table. Stops at the first stop codon, like a real
    ribosome. Any trailing bases that don't form a complete codon are
    ignored.

    Args:
        rna: An RNA string consisting of 'A', 'U', 'C', 'G'.

    Returns:
        The translated protein string, ending at the first stop codon.
    """
    protein = []
    for i in range(0, len(rna) - 2, 3):
        codon = rna[i:i + 3]
        amino_acid = CODON_TABLE.get(codon)
        if amino_acid == 'Stop' or amino_acid is None:
            break
        protein.append(amino_acid)
    return ''.join(protein)


def main():
    import sys
    import os

    input_file = 'prot_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            rna = f.read().strip()
    else:
        rna = sys.stdin.read().strip()

    print(translate(rna))


if __name__ == '__main__':
    main()
