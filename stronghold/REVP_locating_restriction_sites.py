# REVP (Locating Restriction Sites)
# https://rosalind.info/problems/revp/
#
# What a restriction site is:
# Restriction enzymes are bacterial defense proteins that cut DNA at
# specific sequences. Most restriction enzymes are homodimers — two identical
# subunits — and they cut both DNA strands. For the enzyme to cut both
# strands effectively, its target sequence must appear on both strands,
# which happens when the sequence equals its own reverse complement. Such
# sequences are called reverse palindromes (or just palindromes in the
# DNA context).
#
# Why reverse palindromes matter:
# Restriction enzymes recognize these palindromic sequences and cut DNA
# at them. This is the foundation of restriction mapping, cloning, and
# genetic engineering. Different restriction enzymes recognize different
# palindromic sequences (typically 4-8 nucleotides), giving molecular
# biologists a toolkit for cutting DNA at precise locations.
#
# Example: GCATGC is a reverse palindrome because its reverse complement
# is also GCATGC. Reading one strand 5' to 3' gives GCATGC, and reading the
# complementary strand 5' to 3' also gives GCATGC.
#
# What this solution does:
# It scans the DNA string for every substring of length 4 to 12 and checks
# whether it's a reverse palindrome (equals its reverse complement). For
# each match, it reports the 1-based starting position and length.


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


def is_reverse_palindrome(sub):
    """Check if a DNA substring is a reverse palindrome.

    A reverse palindrome equals its reverse complement: each base pairs
    with its complement on the opposite side (A↔T, C↔G).

    Args:
        sub: A DNA substring.

    Returns:
        True if sub is a reverse palindrome, False otherwise.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    n = len(sub)
    for i in range((n + 1) // 2):
        if complement[sub[i]] != sub[n - 1 - i]:
            return False
    return True


def find_reverse_palindromes(dna, min_length=4, max_length=12):
    """Find all reverse palindromes of length 4 to 12 in a DNA string.

    Args:
        dna: A DNA string.
        min_length: Minimum palindrome length (default 4).
        max_length: Maximum palindrome length (default 12).

    Returns:
        A list of (position, length) tuples, where position is 1-based.
    """
    results = []
    for length in range(min_length, max_length + 1):
        for start in range(len(dna) - length + 1):
            sub = dna[start:start + length]
            if is_reverse_palindrome(sub):
                results.append((start + 1, length))  # 1-based position
    return results


def main():
    import sys
    import os

    input_file = 'revp_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            text = f.read()
    else:
        text = sys.stdin.read()

    sequences = parse_fasta(text)
    dna = list(sequences.values())[0]

    palindromes = find_reverse_palindromes(dna)

    for position, length in palindromes:
        print(f'{position} {length}')


if __name__ == '__main__':
    main()
