# ORF (Open Reading Frames)
# https://rosalind.info/problems/orf/
#
# What an open reading frame is:
# An ORF is a stretch of DNA that could potentially code for a protein.
# It starts with a start codon (ATG, which codes for methionine) and ends
# with a stop codon (TAA, TAG, or TGA), with no stop codons in between.
# The ribosome starts translation at the start codon and keeps going until
# it hits a stop codon, at which point it releases the finished protein.
#
# Why six reading frames matter:
# DNA is double-stranded and antiparallel. Either strand can be the coding
# strand, and on each strand there are three ways to group the bases into
# codons (starting at position 0, 1, or 2). That gives six reading frames
# total. The same DNA region can produce completely different proteins
# depending on which frame you read — the triplet grouping shifts everything.
#
# Why ORFs matter in real biology:
# When scientists sequence a genome, they don't know where the genes are.
# One of the first steps in gene finding is to scan for ORFs — long stretches
# without stop codons that start with ATG. Real genes are usually long ORFs,
# while random DNA will hit stop codons frequently. Finding ORFs is how we
# locate candidate genes in newly sequenced DNA.
#
# What this solution does:
# It parses a FASTA DNA string, computes its reverse complement, and scans
# all six reading frames (3 offsets × 2 strands). For each ATG start codon
# found, it translates until the next stop codon and collects the protein
# string. All distinct proteins are printed.


# DNA codon table: maps each DNA codon to its amino acid or 'Stop'.
# Uses T instead of U compared to the RNA codon table.
CODON_TABLE = {
    # Alanine
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    # Arginine
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    # Asparagine
    'AAT': 'N', 'AAC': 'N',
    # Aspartic acid
    'GAT': 'D', 'GAC': 'D',
    # Cysteine
    'TGT': 'C', 'TGC': 'C',
    # Glutamine
    'CAA': 'Q', 'CAG': 'Q',
    # Glutamic acid
    'GAA': 'E', 'GAG': 'E',
    # Glycine
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    # Histidine
    'CAT': 'H', 'CAC': 'H',
    # Isoleucine
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
    # Leucine
    'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    # Lysine
    'AAA': 'K', 'AAG': 'K',
    # Methionine (start codon)
    'ATG': 'M',
    # Phenylalanine
    'TTT': 'F', 'TTC': 'F',
    # Proline
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    # Serine
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
    # Threonine
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    # Tryptophan
    'TGG': 'W',
    # Tyrosine
    'TAT': 'Y', 'TAC': 'Y',
    # Valine
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    # Stop codons
    'TAA': 'Stop', 'TAG': 'Stop', 'TGA': 'Stop',
}

STOP_CODONS = {'TAA', 'TAG', 'TGA'}
START_CODON = 'ATG'


def parse_fasta(text):
    """Parse a FASTA-formatted string into a dictionary of {label: sequence}.

    Args:
        text: A string in FASTA format.

    Returns:
        A dictionary mapping sequence labels (without the '>') to their
        concatenated DNA sequences.
    """
    sequences = {}
    current_label = None
    current_parts = []

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if current_label is not None:
                sequences[current_label] = ''.join(current_parts)
            current_label = line[1:].strip()
            current_parts = []
        else:
            current_parts.append(line)

    if current_label is not None:
        sequences[current_label] = ''.join(current_parts)

    return sequences


def reverse_complement(dna):
    """Return the reverse complement of a DNA string."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(dna))


def translate_orf(dna, start_index):
    """Translate an ORF starting at the given index.

    Reads codons from start_index until a stop codon is reached. Returns
    the protein string (without the stop codon), or None if no stop codon
    is found before the end of the string.

    Args:
        dna: The DNA string.
        start_index: The index of the start codon (ATG).

    Returns:
        The protein string, or None if incomplete.
    """
    protein = []
    i = start_index
    while i + 3 <= len(dna):
        codon = dna[i:i + 3]
        amino_acid = CODON_TABLE.get(codon)
        if amino_acid == 'Stop':
            return ''.join(protein)
        if amino_acid is None:
            return None  # invalid codon, shouldn't happen with valid DNA
        protein.append(amino_acid)
        i += 3
    return None  # no stop codon found


def find_orfs(dna):
    """Find all distinct protein strings from ORFs in all six reading frames.

    Args:
        dna: A DNA string.

    Returns:
        A set of distinct protein strings.
    """
    proteins = set()

    # Check both the forward strand and its reverse complement
    strands = [dna, reverse_complement(dna)]

    for strand in strands:
        # Three reading frames: offset 0, 1, 2
        for offset in range(3):
            i = offset
            while i + 3 <= len(strand):
                codon = strand[i:i + 3]
                if codon == START_CODON:
                    protein = translate_orf(strand, i)
                    if protein is not None and len(protein) > 0:
                        proteins.add(protein)
                i += 3

    return proteins


def main():
    import sys
    import os

    input_file = 'orf_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            text = f.read()
    else:
        text = sys.stdin.read()

    sequences = parse_fasta(text)
    dna = list(sequences.values())[0]  # take the first (and typically only) sequence

    proteins = find_orfs(dna)

    for protein in sorted(proteins):
        print(protein)


if __name__ == '__main__':
    main()
