# CONS (Consensus and Profile)
# https://rosalind.info/problems/cons/
#
# What consensus and profile mean biologically:
# When biologists have several homologous DNA sequences — sequences that
# share a common evolutionary ancestor — they often want to know what the
# ancestral sequence looked like. A consensus string is the best guess at
# that ancestor: at each position, it takes the most common base among all
# the sequences. If most sequences have an 'A' at position 3, the consensus
# has an 'A' there too. The profile matrix shows the full picture: for each
# position, it counts how many times each base (A, C, G, T) appears across
# all the sequences.
#
# Why this matters in real biology:
# Comparing homologous sequences is fundamental to evolutionary biology.
# If you have DNA from several individuals of the same species, or from
# related species, the consensus string tells you what the "typical" sequence
# looks like. Positions where all sequences agree are probably under strong
# selective pressure — they're important enough that mutations there get
# weeded out. Positions with lots of variation are more tolerant of change.
# This is how researchers identify conserved regions, predict functional
# elements, and build phylogenetic trees.
#
# What this solution does:
# It reads a collection of DNA strings in FASTA format, builds the 4×n
# profile matrix (counts of A, C, G, T at each position), derives the
# consensus string by picking the most common base at each position, and
# prints both in the format Rosalind expects.


def parse_fasta(text):
    """Parse a FASTA-formatted string into a dictionary of {label: sequence}.

    FASTA format: each sequence starts with a header line beginning with '>',
    followed by one or more lines of sequence data. This function joins all
    sequence lines for each label into a single string.

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


def consensus_and_profile(sequences):
    """Compute the consensus string and profile matrix from DNA sequences.

    All sequences must be the same length. The profile matrix counts how
    many times each base appears at each position across all sequences.
    The consensus string takes the most common base at each position.

    Args:
        sequences: A list of DNA strings of equal length.

    Returns:
        A tuple of (consensus_string, profile_dict) where profile_dict maps
        each base ('A', 'C', 'G', 'T') to a list of counts per position.
    """
    if not sequences:
        return '', {'A': [], 'C': [], 'G': [], 'T': []}

    length = len(sequences[0])
    bases = ['A', 'C', 'G', 'T']
    profile = {base: [0] * length for base in bases}

    for seq in sequences:
        for i, base in enumerate(seq):
            profile[base][i] += 1

    consensus = []
    for i in range(length):
        best_base = max(bases, key=lambda b: profile[b][i])
        consensus.append(best_base)

    return ''.join(consensus), profile


def main():
    import sys
    import os

    input_file = 'cons_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            text = f.read()
    else:
        text = sys.stdin.read()

    sequences_dict = parse_fasta(text)
    sequences = list(sequences_dict.values())

    consensus, profile = consensus_and_profile(sequences)

    print(consensus)
    for base in ['A', 'C', 'G', 'T']:
        counts = ' '.join(str(c) for c in profile[base])
        print(f'{base}: {counts}')


if __name__ == '__main__':
    main()
