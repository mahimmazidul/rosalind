# ORFR (Finding Genes with ORFs)
# https://rosalind.info/problems/orfr/
#
# What an ORF finder is:
# An ORF finder scans a DNA sequence for open reading frames — regions
# that start with a start codon (ATG) and end at a stop codon (TAA, TAG,
# or TGA) or the end of the sequence, with no stop codons in between.
# It checks all six reading frames (3 on the forward strand, 3 on the
# reverse complement) because either DNA strand can be the coding strand.
#
# Why ORF finding matters:
# ORF finding is one of the first steps in gene annotation. In prokaryotes
# and viruses with compact genomes, ORFs are strong candidates for actual
# genes. Tools like EMBOSS's getorf or SMS 2's ORF finder automate this
# scan, helping researchers identify potential genes in newly sequenced
# DNA.
#
# What this solution does:
# It reads a DNA string, checks all six reading frames for ORFs, translates
# each ORF to protein using Biopython's translate() method, and returns
# the longest protein string found.


from Bio.Seq import Seq
import sys
import os


CODON_TABLE_STANDARD = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

STOP_CODONS = {'TAA', 'TAG', 'TGA'}
START_CODON = 'ATG'


def reverse_complement(dna):
    """Return the reverse complement of a DNA string."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(dna))


def find_orfs_in_strand(dna):
    """Find all ORFs in a single DNA strand.

    Scans all three reading frames. An ORF starts with ATG and ends at
    the first stop codon or at the end of the string.

    Args:
        dna: A DNA string.

    Returns:
        A list of protein strings translated from each ORF.
    """
    proteins = []
    for offset in range(3):
        i = offset
        while i + 3 <= len(dna):
            if dna[i:i + 3] == START_CODON:
                # Translate from this start codon until stop or end
                protein = []
                j = i
                while j + 3 <= len(dna):
                    codon = dna[j:j + 3]
                    aa = CODON_TABLE_STANDARD.get(codon, '*')
                    if aa == '*':
                        break
                    protein.append(aa)
                    j += 3
                if protein:
                    proteins.append(''.join(protein))
                i = j  # skip past this ORF
            else:
                i += 3
    return proteins


def find_longest_orf_protein(dna):
    """Find the longest protein string from any ORF in all 6 reading frames.

    Args:
        dna: A DNA string.

    Returns:
        The longest protein string from any ORF.
    """
    # Check forward strand
    forward_proteins = find_orfs_in_strand(dna)

    # Check reverse complement strand
    rc_dna = reverse_complement(dna)
    reverse_proteins = find_orfs_in_strand(rc_dna)

    all_proteins = forward_proteins + reverse_proteins

    if not all_proteins:
        return ''

    # Return the longest protein
    return max(all_proteins, key=len)


def main():
    input_file = 'orfr_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            lines = f.read().splitlines()
    else:
        lines = sys.stdin.read().splitlines()

    # Filter out empty lines and FASTA headers
    dna_lines = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith('>'):
            continue
        dna_lines.append(line)

    dna = ''.join(dna_lines)

    result = find_longest_orf_protein(dna)
    print(result)


if __name__ == '__main__':
    main()
