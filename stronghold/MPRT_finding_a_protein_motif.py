# MPRT (Finding a Protein Motif)
# https://rosalind.info/problems/mprt/
#
# What a protein motif is:
# A protein motif is a short, conserved sequence of amino acids that plays
# a key role in the protein's structure or function. It might be part of an
# active site where catalysis happens, a binding site where the protein
# attaches to another molecule, or a modification site where a chemical
# group gets added. Motifs tend to stay unchanged through evolution because
# mutating them would break the protein's function.
#
# Why finding motifs matters:
# If you find the same motif in proteins from different species, it's a
# strong sign that the motif has an important function that natural
# selection has preserved. Researchers use motif searches to predict what a
# newly discovered protein does, to find drug targets, and to understand
# how proteins interact with each other.
#
# The N-glycosylation motif:
# N-glycosylation is a common protein modification where a sugar chain
# gets attached to an asparagine (N) residue. The motif N{P}[ST]{P} means:
# asparagine, followed by any amino acid except proline, followed by
# serine or threonine, followed by any amino acid except proline. Proline
# at positions 2 or 4 blocks the modification because proline's rigid ring
# structure changes the protein's shape in a way that prevents the
# glycosylation enzyme from working.
#
# What this solution does:
# It reads UniProt access IDs from a file, fetches each protein's sequence
# from the UniProt database, searches for the N-glycosylation motif, and
# prints the access ID and 1-based positions for each protein that has at
# least one match.


def fetch_protein_sequence(access_id):
    """Fetch a protein's amino acid sequence from UniProt in FASTA format.

    Args:
        access_id: A UniProt access ID (e.g., 'B5ZC00' or
                   'P07204_TRBM_HUMAN').

    Returns:
        The protein sequence as a string, or None if the fetch fails.
    """
    import urllib.request
    import urllib.error

    # Strip any suffix after underscore (e.g., P07204_TRBM_HUMAN -> P07204)
    accession = access_id.split('_')[0]
    url = f'https://rest.uniprot.org/uniprotkb/{accession}.fasta'
    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            fasta_text = response.read().decode('utf-8')
    except (urllib.error.URLError, Exception):
        return None

    # Parse FASTA: skip header lines, join sequence lines
    sequence_parts = []
    for line in fasta_text.splitlines():
        if line.startswith('>') or not line.strip():
            continue
        sequence_parts.append(line.strip())

    return ''.join(sequence_parts)


def find_n_glycosylation_motifs(protein):
    """Find all 1-based positions of the N-glycosylation motif N{P}[ST]{P}.

    The motif is: N, then any amino acid except P, then S or T,
    then any amino acid except P.

    Args:
        protein: A string of amino acid single-letter codes.

    Returns:
        A list of 1-based starting positions where the motif occurs.
    """
    positions = []
    for i in range(len(protein) - 3):
        if (
            protein[i] == 'N'
            and protein[i + 1] != 'P'
            and protein[i + 2] in ('S', 'T')
            and protein[i + 3] != 'P'
        ):
            positions.append(i + 1)  # 1-based position
    return positions


def main():
    import sys
    import os

    input_file = 'mprt_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            access_ids = [line.strip() for line in f if line.strip()]
    else:
        access_ids = [line.strip() for line in sys.stdin if line.strip()]

    for access_id in access_ids:
        protein = fetch_protein_sequence(access_id)
        if protein is None:
            continue
        positions = find_n_glycosylation_motifs(protein)
        if positions:
            print(access_id)
            print(' '.join(str(p) for p in positions))


if __name__ == '__main__':
    main()
