# SUBO (Suboptimal Local Alignment)
# https://rosalind.info/problems/subo/
#
# What suboptimal local alignment is:
# Standard local alignment (like Smith-Waterman) finds only the single best
# matching region between two sequences. But biological sequences often
# contain multiple similar regions — for example, transposable elements
# ("jumping genes") that have copied themselves to multiple locations in
# the genome. These copies are usually not identical; they accumulate
# mutations over time. Suboptimal alignment finds multiple similar regions
# simultaneously, not just the best one.
#
# Why this matters in real biology:
# Transposons and other repetitive elements make up a large fraction of
# many genomes. They can cause gene duplications, deletions, and
# rearrangements. Finding all copies of a repeat — even imperfect ones —
# is essential for understanding genome structure and evolution. Tools
# like Lalign (from the FASTA package) implement suboptimal alignment to
# find these multiple similar regions.
#
# What this solution does:
# Since EMBOSS Lalign is not available locally, this solution approximates
# its function. It finds a short repeat (32-40 bp) shared by two DNA
# strings, allowing up to 3 differences (substitutions or indels) per
# occurrence, then counts how many times the repeat appears in each string.
# It uses a k-mer seed approach: extract candidate repeats from the first
# string, then count approximate occurrences of each candidate in both
# strings using edit distance.


def edit_distance(s1, s2):
    """Compute the edit distance (Levenshtein) between two strings.

    Uses dynamic programming with O(min(len(s1), len(s2))) space.

    Args:
        s1: First string.
        s2: Second string.

    Returns:
        The edit distance (minimum number of insertions, deletions,
        and substitutions to transform s1 into s2).
    """
    if len(s1) < len(s2):
        s1, s2 = s2, s1

    prev = list(range(len(s2) + 1))
    for i in range(1, len(s1) + 1):
        curr = [i] + [0] * len(s2)
        for j in range(1, len(s2) + 1):
            cost = 0 if s1[i - 1] == s2[j - 1] else 1
            curr[j] = min(
                prev[j] + 1,      # deletion
                curr[j - 1] + 1,  # insertion
                prev[j - 1] + cost  # substitution
            )
        prev = curr
    return prev[-1]


def count_approximate_occurrences(text, pattern, max_diff=3):
    """Count non-overlapping occurrences of pattern in text with at most
    max_diff differences.

    Slides the pattern across the text. When a match is found, skips ahead
    by the pattern length to avoid counting overlapping occurrences.

    Args:
        text: The text to search within.
        pattern: The pattern to search for.
        max_diff: Maximum allowed edit distance (default 3).

    Returns:
        The number of non-overlapping approximate occurrences.
    """
    count = 0
    plen = len(pattern)
    i = 0
    while i <= len(text) - plen:
        substring = text[i:i + plen]
        if edit_distance(substring, pattern) <= max_diff:
            count += 1
            i += plen  # skip past this match (non-overlapping)
        else:
            i += 1
    return count


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


def find_repeat(s1, s2, min_len=32, max_len=40, max_diff=3):
    """Find a shared inexact repeat between two DNA strings.

    Args:
        s1: First DNA string.
        s2: Second DNA string.
        min_len: Minimum repeat length.
        max_len: Maximum repeat length.
        max_diff: Maximum differences per occurrence.

    Returns:
        A tuple of (count_in_s1, count_in_s2) or None if no repeat found.
    """
    # Try different repeat lengths
    for length in range(min_len, max_len + 1):
        best_pattern = None
        best_total = 0

        # Extract all k-mers of this length from s1
        for i in range(len(s1) - length + 1):
            candidate = s1[i:i + length]

            # Count approximate occurrences in both strings
            count1 = count_approximate_occurrences(s1, candidate, max_diff)
            count2 = count_approximate_occurrences(s2, candidate, max_diff)

            # The repeat must appear at least twice in each string
            if count1 >= 2 and count2 >= 2:
                total = count1 + count2
                if total > best_total:
                    best_total = total
                    best_pattern = candidate

        if best_pattern is not None:
            count1 = count_approximate_occurrences(s1, best_pattern, max_diff)
            count2 = count_approximate_occurrences(s2, best_pattern, max_diff)
            return count1, count2

    return None


def main():
    import sys
    import os

    input_file = 'subo_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            text = f.read()
    else:
        text = sys.stdin.read()

    sequences = parse_fasta(text)
    s1 = sequences[0][1]
    s2 = sequences[1][1]

    result = find_repeat(s1, s2)
    if result:
        print(f'{result[0]} {result[1]}')
    else:
        print("No repeat found")


if __name__ == '__main__':
    main()
