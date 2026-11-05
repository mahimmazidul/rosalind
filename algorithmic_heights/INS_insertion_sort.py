# INS (Insertion Sort)
# https://rosalind.info/problems/ins/
#
# What the algorithm is:
# Insertion sort builds the final sorted array one element at a time. It
# takes each element from the input and inserts it into its correct position
# in the already-sorted portion of the array, shifting larger elements one
# position to the right to make room.
#
# How it works conceptually:
# Think of sorting a hand of playing cards. You pick up one card at a time
# and insert it into the correct position among the cards you've already
# sorted. The cards to the right of your insertion point shift over to make
# room.
#
# Time and space complexity:
# Time: O(n^2) worst case — each element may need to be compared and shifted
# against all previously sorted elements.
# Space: O(1) — sorts in place, using only a few integer variables.
#
# Where this appears in bioinformatics pipelines:
# Insertion sort is efficient for small or nearly sorted datasets. It's
# often used as the base case in divide-and-conquer sorts (like merge sort)
# when the subarray size drops below a threshold. In bioinformatics, it can
# be useful for sorting small lists of variants or annotations that are
# already mostly in order.


def insertion_sort_swaps(arr):
    """Count the number of swaps performed by insertion sort.

    Each time an element is shifted one position to the right to make
    room for the key being inserted, that counts as one swap.

    Args:
        arr: A list of integers.

    Returns:
        The total number of swaps performed.
    """
    swaps = 0
    for i in range(1, len(arr)):
        key = arr[i]
        j = i - 1
        while j >= 0 and arr[j] > key:
            arr[j + 1] = arr[j]
            swaps += 1
            j -= 1
        arr[j + 1] = key
    return swaps


def main():
    import sys
    import os

    input_file = 'ins_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            content = f.read()
    else:
        content = sys.stdin.read()

    lines = content.strip().splitlines()
    n = int(lines[0].strip())
    arr = list(map(int, lines[1].strip().split()))

    swaps = insertion_sort_swaps(arr)
    print(swaps)


if __name__ == '__main__':
    main()
