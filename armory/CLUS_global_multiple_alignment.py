# CLUS (Global Multiple Alignment)
# https://rosalind.info/problems/clus/
#
# What multiple sequence alignment (MSA) is:
# MSA aligns three or more biological sequences simultaneously, inserting
# gaps so that homologous positions line up in columns. It's used to
# identify conserved regions, infer evolutionary relationships, and detect
# outliers — sequences that don't fit the pattern of the others.
#
# Why MSA matters in real bioinformatics:
# MSA is fundamental to comparative genomics. It reveals conserved domains,
# helps build phylogenetic trees, and can identify sequences that are
# contaminants or mislabeled samples — for example, detecting tilapia
# passed off as salmon in food safety testing.
#
# The Clustal approach:
# ClustalW uses a progressive alignment strategy: (1) compute all pairwise
# alignments to build a distance matrix, (2) build a guide tree from the
# matrix, (3) align sequences progressively following the tree. The most
# different sequence is the one with the greatest average distance from
# all others.
#
# What this solution does:
# Since ClustalW is not available locally, this solution approximates it
# using Biopython's PairwiseAligner. It computes pairwise global alignment
# distances between all pairs of DNA strings, then identifies the string
# with the highest average distance to all others — the most different
# from the group.


from Bio.Align import PairwiseAligner, substitution_matrices
from Bio import SeqIO
import sys
import os


def parse_fasta(text):
    """Parse a FASTA-formatted string into a list of (label, sequence) tuples.

    Args:
        text: A string in FASTA format.

    Returns:
        A list of (label, sequence) tuples.
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


def alignment_distance(seq1, seq2, aligner):
    """Compute the distance between two sequences from their global alignment.

    Distance = 1 - (matches / alignment_length), where matches are positions
    where both sequences have the same non-gap residue.

    Args:
        seq1: First DNA string.
        seq2: Second DNA string.
        aligner: A configured PairwiseAligner.

    Returns:
        A distance value between 0 (identical) and 1 (completely different).
    """
    alignments = aligner.align(seq1, seq2)
    aln = alignments[0]

    aligned1 = str(aln[0])
    aligned2 = str(aln[1])

    matches = 0
    aligned_length = 0
    for a, b in zip(aligned1, aligned2):
        if a == '-' or b == '-':
            continue  # skip gap positions for match counting
        aligned_length += 1
        if a == b:
            matches += 1

    if aligned_length == 0:
        return 1.0

    return 1.0 - (matches / aligned_length)


def find_most_different(sequences, aligner):
    """Find the sequence most different from all others.

    Computes pairwise distances between all pairs, then for each sequence
    computes its average distance to all others. Returns the label of the
    sequence with the highest average distance.

    Args:
        sequences: A list of (label, sequence) tuples.
        aligner: A configured PairwiseAligner.

    Returns:
        The label of the most different sequence.
    """
    n = len(sequences)
    if n <= 1:
        return sequences[0][0] if sequences else None

    # Compute average distance for each sequence
    avg_distances = []
    for i in range(n):
        total_distance = 0.0
        count = 0
        for j in range(n):
            if i == j:
                continue
            dist = alignment_distance(sequences[i][1], sequences[j][1], aligner)
            total_distance += dist
            count += 1
        avg_distances.append((sequences[i][0], total_distance / count))

    # Return the label with the highest average distance
    avg_distances.sort(key=lambda x: x[1], reverse=True)
    return avg_distances[0][0]


def main():
    input_file = 'clus_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            text = f.read()
    else:
        text = sys.stdin.read()

    sequences = parse_fasta(text)

    # Configure the aligner for global alignment with DNAfull-like scoring
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load('NUC.4.4')
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1
    aligner.mode = 'global'

    result = find_most_different(sequences, aligner)
    print(result)


if __name__ == '__main__':
    main()
