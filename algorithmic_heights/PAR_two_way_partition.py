# PAR (Two Way Partition)
# https://rosalind.info/problems/par/
#
# What the algorithm is:
# Two-way partition is the core subroutine of Quick Sort. It rearranges an
# array so that all elements less than or equal to the pivot (the first
# element) come before it, and all elements greater than the pivot come
# after it. The pivot ends up in its final sorted position.
#
# How it works conceptually:
# Use two pointers scanning from both ends. The left pointer finds elements
# greater than the pivot, the right pointer finds elements less than or
# equal to the pivot. When both are found, swap them. Repeat until the
# pointers cross, then place the pivot in its correct position.
#
# Time and space complexity:
# Time: O(n) — single pass through the array.
# Space: O(1) — in-place, only a few integer variables.
#
# Where this appears in bioinformatics pipelines:
# Partitioning is used in Quick Sort for sorting large genomic datasets.
# It's also the basis for selection algorithms (finding the kth smallest
# element), useful for finding median values in genomics statistics.


def two_way_partition(arr):
    """Partition arr around the first element (pivot).

    Elements <= pivot go to the left, pivot goes in the middle,
    elements > pivot go to the right.

    Args:
        arr: A list of integers.

    Returns:
        The permuted list with the pivot in its correct position.
    """
    if len(arr) <= 1:
        return arr

    pivot = arr[0]
    left = []
    right = []

    for x in arr[1:]:
        if x <= pivot:
            left.append(x)
        else:
            right.append(x)

    return left + [pivot] + right


def main():
    import sys
    import os

    input_file = 'par_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            content = f.read()
    else:
        content = sys.stdin.read()

    lines = content.strip().splitlines()
    n = int(lines[0].strip())
    arr = list(map(int, lines[1].strip().split()))

    result = two_way_partition(arr)
    print(' '.join(str(x) for x in result))


if __name__ == '__main__':
    main()
