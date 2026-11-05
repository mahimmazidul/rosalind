# PRTM (Calculating Protein Mass)
# https://rosalind.info/problems/prtm/
#
# What protein mass is:
# Every amino acid has a characteristic mass. When amino acids link together
# to form a protein, they form peptide bonds and lose water molecules. The
# mass of a protein is the sum of the masses of its amino acid residues
# (amino acids minus water) plus the mass of one water molecule at the ends.
#
# Why protein mass matters:
# Mass spectrometry is a key tool in proteomics. It measures the mass of
# proteins and peptides with high precision. By comparing measured masses
# to calculated masses, researchers can identify proteins, detect
# post-translational modifications, and discover mutations. This is how
# mass spectrometry-based protein identification works.
#
# Monoisotopic vs average mass:
# Each element has multiple isotopes with slightly different masses.
# Monoisotopic mass uses the most abundant isotope of each element and is
# the standard in mass spectrometry. Average mass weighs isotopes by their
# natural abundance. For precise mass spectrometry work, monoisotopic mass
# is used.
#
# What this solution does:
# It looks up the monoisotopic mass of each amino acid in the protein
# string and sums them up. The result is the total monoisotopic mass of
# the protein.


# Monoisotopic masses of the 20 standard amino acids (in Daltons).
# Source: standard monoisotopic mass table used in mass spectrometry.
AMINO_ACID_MASSES = {
    'A': 71.03711,   # Alanine
    'C': 103.00919,  # Cysteine
    'D': 115.02694,  # Aspartic acid
    'E': 129.04259,  # Glutamic acid
    'F': 147.06841,  # Phenylalanine
    'G': 57.02146,   # Glycine
    'H': 137.05891,  # Histidine
    'I': 113.08406,  # Isoleucine
    'K': 128.09496,  # Lysine
    'L': 113.08406,  # Leucine
    'M': 131.04049,  # Methionine
    'N': 114.04293,  # Asparagine
    'P': 97.05276,   # Proline
    'Q': 128.05858,  # Glutamine
    'R': 156.10111,  # Arginine
    'S': 87.03203,   # Serine
    'T': 101.04768,  # Threonine
    'V': 99.06841,   # Valine
    'W': 186.07931,  # Tryptophan
    'Y': 163.06333,  # Tyrosine
}


def calculate_protein_mass(protein):
    """Calculate the total monoisotopic mass of a protein string.

    Args:
        protein: A string of amino acid single-letter codes.

    Returns:
        The total monoisotopic mass as a float.
    """
    return sum(AMINO_ACID_MASSES[aa] for aa in protein)


def main():
    import sys
    import os

    input_file = 'prtm_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            protein = f.read().strip()
    else:
        protein = sys.stdin.read().strip()

    mass = calculate_protein_mass(protein)
    print(f"{mass:.3f}")


if __name__ == '__main__':
    main()
