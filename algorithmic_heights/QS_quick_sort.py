# QS (Quick Sort)
# https://rosalind.info/problems/qs/
#
# What the algorithm is:
# Quicksort — a divide-and-conquer sorting algorithm that partitions the
# array around a pivot (elements smaller than pivot go left, larger go
# right), then recursively sorts the two subarrays.
#
# How it works conceptually:
# Pick a pivot, partition the array so everything less than the pivot is
# on its left and everything greater is on its right. The pivot is now in
# its final sorted position. Recurse on the left and right parts.
# This implementation is iterative (explicit stack) to avoid Python's
# recursion limit on large inputs (n up to 10^5).
#
# Time and space complexity:
# Time: O(n log n) average, O(n^2) worst case.
# Space: O(log n) average stack depth with random pivots.
#
# Graph representation: N/A (array sorting).
#
# Where this appears in bioinformatics:
# Quicksort underlies sorting of large genomic datasets — e.g. sorting
# aligned reads by position, sorting variant calls by coordinate, and
# any pipeline step that orders millions of short sequences.

import sys
import os
import random


def quicksort(arr):
    """Sort arr in place using iterative quicksort with random pivots."""
    if len(arr) <= 1:
        return arr

    # Explicit stack of (low, high) subarray bounds
    stack = [(0, len(arr) - 1)]

    while stack:
        low, high = stack.pop()
        if low >= high:
            continue

        # Random pivot to avoid worst-case on sorted input
        pivot_idx = random.randint(low, high)
        arr[pivot_idx], arr[high] = arr[high], arr[pivot_idx]
        pivot = arr[high]

        # Partition: elements <= pivot go to the left
        i = low
        for j in range(low, high):
            if arr[j] <= pivot:
                arr[i], arr[j] = arr[j], arr[i]
                i += 1
        arr[i], arr[high] = arr[high], arr[i]

        # Push larger subarray first so the stack stays shallow
        left_size = i - 1 - low
        right_size = high - (i + 1)
        if left_size > right_size:
            stack.append((low, i - 1))
            stack.append((i + 1, high))
        else:
            stack.append((i + 1, high))
            stack.append((low, i - 1))

    return arr


def main():
    input_file = 'qs_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            content = f.read()
    else:
        content = sys.stdin.read()

    lines = content.strip().splitlines()
    n = int(lines[0].strip())
    arr = list(map(int, lines[1].split()))

    sorted_arr = quicksort(arr)
    print(' '.join(str(x) for x in sorted_arr))


if __name__ == '__main__':
    main()
