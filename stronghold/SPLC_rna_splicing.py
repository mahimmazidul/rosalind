# SPLC (RNA Splicing)
# https://rosalind.info/problems/splc/
#
# What RNA splicing is:
# In eukaryotic cells, genes are not continuous stretches of DNA. They are
# interrupted by non-coding sequences called introns. When a gene is
# transcribed, the entire region (introns and exons) is copied into a
# pre-mRNA. Then the spliceosome cuts out the introns and glues the exons
# together to form mature mRNA, which is then translated into protein.
#
# Why splicing matters:
# Splicing allows one gene to code for multiple proteins through
# alternative splicing — different combinations of exons can be included
# or excluded. It also means that mutations in intron-exon boundaries can
# cause disease by disrupting the splicing process. Understanding splicing
# is essential for interpreting genomes and designing gene therapies.
#
# What this solution does:
# It reads a FASTA file where the first sequence is the genomic DNA and the
# remaining sequences are introns. It removes each intron from the DNA,
# transcribes the result to RNA (T→U), translates using the standard
# genetic code (stopping at the first stop codon), and prints the protein.


# Standard genetic code: maps each RNA codon to its amino acid or 'Stop'.
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


def parse_fasta(text):
    """Parse a FASTA-formatted string into a list of (label, sequence) tuples.

    Args:
        text: A string in FASTA format.

    Returns:
        A list of (label, sequence) tuples in the order they appear.
    """
    sequences = []
    current_label = None
    current_parts = []

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if current_label is not None:
                sequences.append((current_label, ''.join(current_parts)))
            current_label = line[1:].strip()
            current_parts = []
        else:
            current_parts.append(line)

    if current_label is not None:
        sequences.append((current_label, ''.join(current_parts)))

    return sequences


def remove_introns(dna, introns):
    """Remove all intron sequences from the DNA string.

    Args:
        dna: The genomic DNA string.
        introns: A list of intron strings to remove.

    Returns:
        The DNA string with all introns removed.
    """
    for intron in introns:
        dna = dna.replace(intron, '')
    return dna


def transcribe(dna):
    """Transcribe DNA to RNA by replacing T with U."""
    return dna.replace('T', 'U')


def translate(rna):
    """Translate an RNA string into a protein string.

    Reads codons from the start until a stop codon is reached.

    Args:
        rna: An RNA string.

    Returns:
        The protein string.
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

    input_file = 'splc_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            text = f.read()
    else:
        text = sys.stdin.read()

    sequences = parse_fasta(text)

    # First sequence is the genomic DNA, the rest are introns
    dna = sequences[0][1]
    introns = [seq for _, seq in sequences[1:]]

    # Remove introns, transcribe, translate
    exons_dna = remove_introns(dna, introns)
    rna = transcribe(exons_dna)
    protein = translate(rna)

    print(protein)


if __name__ == '__main__':
    main()
