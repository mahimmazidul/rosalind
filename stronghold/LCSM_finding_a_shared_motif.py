# LCSM (Finding a Shared Motif)
# https://rosalind.info/problems/lcsm/
#
# What a shared motif is:
# When biologists compare DNA sequences from different organisms, they look
# for regions that are identical across all of them. A shared motif is a
# substring that appears in every sequence in the collection. The longer
# the shared motif, the more likely it represents something biologically
# important — a gene, a regulatory region, or a functional RNA element that
# has been conserved through evolution because it matters.
#
# Why finding the longest shared motif matters:
# If two species share a long stretch of identical DNA, it's very unlikely
# to be a coincidence. It suggests they inherited that sequence from a
# common ancestor and that natural selection has preserved it because
# changing it would be harmful. Finding these shared regions helps
# researchers identify genes, understand evolutionary relationships, and
# pinpoint functionally important parts of the genome.
#
# Why the naive approach works here:
# The simplest approach is to take the shortest string and try all of its
# substrings from longest to shortest, checking whether each one appears in
# every other string. With at most 100 strings of at most 1000 bases each,
# this brute-force approach is fast enough. More sophisticated algorithms
# (like suffix trees) exist for much larger inputs, but for Rosalind's
# constraints, the straightforward method gets the job done.


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


def longest_common_substring(sequences):
    """Find a longest common substring shared by all sequences.

    Takes the shortest sequence and generates all its substrings from
    longest to shortest, checking whether each appears in every other
    sequence. Returns the first match found.

    Args:
        sequences: A dictionary of {label: sequence}.

    Returns:
        A string that is a longest common substring of all sequences.
    """
    seqs = list(sequences.values())
    if not seqs:
        return ''

    shortest = min(seqs, key=len)
    others = [s for s in seqs if s is not shortest]

    # Try substrings from longest to shortest
    for length in range(len(shortest), 0, -1):
        for start in range(len(shortest) - length + 1):
            candidate = shortest[start:start + length]
            if all(candidate in other for other in others):
                return candidate

    return ''


def main():
    import sys
    import os

    input_file = 'lcsm_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            text = f.read()
    else:
        text = sys.stdin.read()

    sequences = parse_fasta(text)
    result = longest_common_substring(sequences)
    print(result)


if __name__ == '__main__':
    main()
