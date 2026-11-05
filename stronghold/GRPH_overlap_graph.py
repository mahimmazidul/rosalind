# GRPH (Overlap Graphs)
# https://rosalind.info/problems/grph/
#
# What an overlap graph is:
# An overlap graph is a directed graph where each node represents a DNA
# string, and there's a directed edge from string s to string t when a
# suffix of s matches a prefix of t. In this problem, we use k=3, so we look
# for length-3 suffixes matching length-3 prefixes. We exclude self-loops
# (a string doesn't connect to itself).
#
# Why overlap graphs matter in real biology:
# Overlap graphs are fundamental to genome assembly. When you sequence a
# genome, you don't get one long continuous string — you get many short
# fragments called reads. To reconstruct the original genome, you need to
# figure out how these reads overlap with each other. If the end of read A
# matches the beginning of read B, they probably come from adjacent parts of
# the genome. Building an overlap graph is the first step in many genome
# assembly algorithms: each read is a node, and edges represent potential
# overlaps. A path through this graph represents a possible reconstruction
# of the genome.
#
# What this solution does:
# It reads DNA strings in FASTA format, then for every ordered pair of
# distinct strings (s, t), checks whether the last 3 characters of s equal
# the first 3 characters of t. When they match, it prints the pair as a
# directed edge in the adjacency list.


def parse_fasta(text):
    """Parse a FASTA-formatted string into a dictionary of {label: sequence}.

    FASTA format: each sequence starts with a header line beginning with '>',
    followed by one or more lines of sequence data.

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


def overlap_graph(sequences, k=3):
    """Build the overlap graph O_k for a collection of DNA strings.

    For every ordered pair of distinct strings (s, t), adds a directed edge
    s -> t when the length-k suffix of s matches the length-k prefix of t.

    Args:
        sequences: A dictionary of {label: sequence}.
        k: The overlap length (default 3).

    Returns:
        A list of (source_label, target_label) tuples representing edges.
    """
    edges = []
    labels = list(sequences.keys())
    for i, s_label in enumerate(labels):
        s = sequences[s_label]
        for j, t_label in enumerate(labels):
            if i == j:
                continue  # no self-loops
            t = sequences[t_label]
            if len(s) >= k and len(t) >= k and s[-k:] == t[:k]:
                edges.append((s_label, t_label))
    return edges


def main():
    import sys
    import os

    input_file = 'grph_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            text = f.read()
    else:
        text = sys.stdin.read()

    sequences = parse_fasta(text)
    edges = overlap_graph(sequences, k=3)

    for source, target in edges:
        print(f'{source} {target}')


if __name__ == '__main__':
    main()
