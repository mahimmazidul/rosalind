# NEED (Pairwise Global Alignment)
# https://rosalind.info/problems/need/
#
# What EMBOSS Needle is:
# EMBOSS Needle is a tool for performing global pairwise sequence alignment.
# It takes two sequences (DNA, RNA, or protein) and finds the optimal global
# alignment — meaning it aligns the entire length of both sequences, inserting
# gaps where needed to maximize the alignment score. Needle uses the
# Needleman-Wunsch algorithm, a dynamic programming approach that guarantees
# finding the optimal alignment given a scoring scheme.
#
# Why global alignment matters in real bioinformatics:
# Global alignment is used when two sequences are expected to be similar
# across their entire length — for example, two homologous genes from
# different species, or two reads from the same genomic region. The alignment
# score quantifies how similar they are, which helps in identifying
# homologous sequences, inferring evolutionary relationships, and detecting
# conserved functional regions.
#
# The scoring scheme:
# This solution uses the DNAfull scoring matrix (match = +5, mismatch = -4),
# which is equivalent to Biopython's NUC.4.4 matrix. Gap opening penalty is
# 10 and gap extension penalty is 1 (affine gap penalty). Affine gap
# penalties reflect the biological reality that a single multi-base insertion
# or deletion is more likely than several independent single-base events.
#
# What this solution does:
# Since EMBOSS is not available locally, this solution uses Biopython's
# PairwiseAligner with the NUC.4.4 substitution matrix as a substitute for
# EMBOSS Needle's DNAfull matrix. It fetches DNA sequences from GenBank using
# Entrez, performs a global alignment, and returns the maximum alignment
# score.


from Bio import Entrez
from Bio.Align import PairwiseAligner, substitution_matrices

Entrez.email = 'mahimmazidul@users.noreply.github.com'


def fetch_dna_sequence(access_id):
    """Fetch a DNA sequence from GenBank by accession ID.

    Args:
        access_id: A GenBank accession ID (e.g., 'JX205496.1').

    Returns:
        The DNA sequence as a string, or None if the fetch fails.
    """
    try:
        handle = Entrez.efetch(db='nucleotide', id=access_id, rettype='fasta')
        fasta_text = handle.read()
        handle.close()
        lines = fasta_text.strip().splitlines()
        return ''.join(line.strip() for line in lines if not line.startswith('>'))
    except Exception:
        return None


def global_alignment_score(seq1, seq2):
    """Compute the maximum global alignment score between two DNA strings.

    Uses the NUC.4.4 substitution matrix (equivalent to EMBOSS DNAfull:
    match = +5, mismatch = -4) with affine gap penalties
    (open = 10, extend = 1).

    Args:
        seq1: First DNA string.
        seq2: Second DNA string.

    Returns:
        The maximum global alignment score as an integer.
    """
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load('NUC.4.4')
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1
    aligner.mode = 'global'

    alignments = aligner.align(seq1, seq2)
    return int(alignments[0].score)


def main():
    import sys
    import os

    input_file = 'need_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            access_ids = f.read().split()
    else:
        access_ids = sys.stdin.read().split()

    if len(access_ids) < 2:
        print("Need two GenBank IDs")
        return

    seq1 = fetch_dna_sequence(access_ids[0])
    seq2 = fetch_dna_sequence(access_ids[1])

    if seq1 is None or seq2 is None:
        print("Failed to fetch sequences")
        return

    score = global_alignment_score(seq1, seq2)
    print(score)


if __name__ == '__main__':
    main()
