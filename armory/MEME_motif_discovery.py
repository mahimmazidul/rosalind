# MEME (Motif Discovery)
# https://rosalind.info/problems/meme/
#
# What MEME is:
# MEME (Multiple EM for Motif Elicitation) is a tool for discovering novel
# motifs in a collection of DNA or protein sequences. A motif is a short,
# conserved sequence pattern that likely has a biological function — such as
# a protein binding site, a modification site, or a structural element.
# MEME represents each motif as a position-specific probability matrix and
# assigns it an E-value indicating statistical significance.
#
# Why MEME exists in real bioinformatics workflows:
# Databases like ProSite only catalog motifs that have already been
# discovered and characterized. When you have a new set of sequences and
# suspect they share a regulatory or functional pattern, you need a de novo
# motif discovery tool. MEME finds these patterns without prior knowledge,
# which is essential for studying newly discovered genes or non-coding
# regulatory regions.
#
# What this solution does:
# Since the full MEME suite is an external tool not available in this
# environment, this solution approximates MEME's output by finding a shared
# motif region across all input sequences. It locates a conserved anchor
# k-mer common to all sequences, extends it left and right while all
# sequences have residues at each position, collects the residue set at
# each position, and outputs a regular expression (using [XY] notation for
# variable positions, as MEME does for motif regular expressions).


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


def find_common_kmer(sequences, k=10):
    """Find a k-mer of length k that appears in all sequences.

    Uses the shortest sequence as the source of candidate kmers.

    Args:
        sequences: A list of (label, sequence) tuples.
        k: Length of the k-mer to search for.

    Returns:
        A tuple of (kmer, positions) where positions is a list of the
        starting index of the kmer in each sequence, or (None, None) if
        no common kmer is found.
    """
    seqs = [seq for _, seq in sequences]
    shortest = min(seqs, key=len)
    shortest_idx = seqs.index(shortest)

    for i in range(len(shortest) - k + 1):
        kmer = shortest[i:i + k]
        positions = [-1] * len(seqs)
        positions[shortest_idx] = i
        found_all = True
        for j, seq in enumerate(seqs):
            if j == shortest_idx:
                continue
            pos = seq.find(kmer)
            if pos == -1:
                found_all = False
                break
            positions[j] = pos
        if found_all:
            return kmer, positions

    return None, None


def extend_motif(sequences, anchor, anchor_positions):
    """Extend the anchor motif left and right while all sequences agree.

    At each position, collects the set of residues found across all
    sequences at the aligned position. Single residues are output as-is;
    multiple residues are output as [XYZ] character classes.

    Args:
        sequences: A list of (label, sequence) tuples.
        anchor: The conserved anchor k-mer string.
        anchor_positions: List of the anchor's start index in each sequence.

    Returns:
        A regex-style motif string.
    """
    seqs = [seq for _, seq in sequences]
    anchor_len = len(anchor)
    n_seqs = len(seqs)

    # Extend left, but only while all sequences share the same residue
    # (or at most 2 distinct residues, which becomes a character class)
    left_ext = 0
    while True:
        left_ext += 1
        residues = set()
        all_have = True
        for i in range(n_seqs):
            pos = anchor_positions[i] - left_ext
            if pos < 0:
                all_have = False
                break
            residues.add(seqs[i][pos])
        if not all_have or len(residues) > 2:
            left_ext -= 1
            break

    # Extend right with the same conservation constraint
    right_ext = 0
    while True:
        right_ext += 1
        residues = set()
        all_have = True
        for i in range(n_seqs):
            pos = anchor_positions[i] + anchor_len + right_ext - 1
            if pos >= len(seqs[i]):
                all_have = False
                break
            residues.add(seqs[i][pos])
        if not all_have or len(residues) > 2:
            right_ext -= 1
            break

    # Build the motif pattern
    start_positions = [anchor_positions[i] - left_ext for i in range(n_seqs)]
    motif_len = left_ext + anchor_len + right_ext

    pattern_parts = []
    for j in range(motif_len):
        residues = []
        seen = set()
        for i in range(n_seqs):
            pos = start_positions[i] + j
            r = seqs[i][pos]
            if r not in seen:
                seen.add(r)
                residues.append(r)
        if len(residues) == 1:
            pattern_parts.append(residues[0])
        else:
            pattern_parts.append('[' + ''.join(residues) + ']')

    return ''.join(pattern_parts)


def main():
    import sys
    import os

    input_file = 'meme_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            text = f.read()
    else:
        text = sys.stdin.read()

    sequences = parse_fasta(text)

    # Try progressively shorter kmers until we find a common one
    anchor = None
    anchor_positions = None
    for k in range(15, 5, -1):
        anchor, anchor_positions = find_common_kmer(sequences, k)
        if anchor is not None:
            break

    if anchor is None:
        # No common kmer found; output nothing
        return

    motif = extend_motif(sequences, anchor, anchor_positions)
    print(motif)


if __name__ == '__main__':
    main()
