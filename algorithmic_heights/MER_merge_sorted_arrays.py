# MER (Merge Two Sorted Arrays)
# https://rosalind.info/problems/mer/
#
# What the algorithm is:
# Merging combines two sorted arrays into a single sorted array. It's the
# core subroutine of merge sort. The algorithm uses two pointers — one for
# each array — and repeatedly picks the smaller of the two current elements,
# appending it to the result and advancing that pointer.
#
# How it works conceptually:
# The smallest element of the merged array must be either the first element
# of A or the first element of B, whichever is smaller. Pick it, then repeat
# with the remaining elements. This is like merging two sorted piles of
# cards: look at the top card of each pile, take the smaller one, and
# continue.
#
# Time and space complexity:
# Time: O(n + m) — each element is examined exactly once.
# Space: O(n + m) — the output array contains all n + m elements.
#
# Where this appears in bioinformatics pipelines:
# Merging sorted lists is fundamental when combining results from multiple
# sources — for example, merging sorted lists of variants from different
# chromosomes, or combining sorted alignment positions from multiple
# sequencing runs. External sort algorithms for large genomic datasets also
# rely heavily on k-way merging.


def merge_sorted_arrays(A, B):
    """Merge two sorted arrays into a single sorted array.

    Args:
        A: A sorted list of integers.
        B: A sorted list of integers.

    Returns:
        A sorted list containing all elements of A and B.
    """
    i, j = 0, 0
    merged = []

    while i < len(A) and j < len(B):
        if A[i] <= B[j]:
            merged.append(A[i])
            i += 1
        else:
            merged.append(B[j])
            j += 1

    # Append any remaining elements
    merged.extend(A[i:])
    merged.extend(B[j:])

    return merged


def main():
    import sys
    import os

    input_file = 'mer_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            content = f.read()
    else:
        content = sys.stdin.read()

    lines = content.strip().splitlines()
    n = int(lines[0].strip())
    A = list(map(int, lines[1].strip().split()))
    m = int(lines[2].strip())
    B = list(map(int, lines[3].strip().split()))

    merged = merge_sorted_arrays(A, B)
    print(' '.join(str(x) for x in merged))


if __name__ == '__main__':
    main()
