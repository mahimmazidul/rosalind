# MED (Median / Selection)
# https://rosalind.info/problems/med/
#
# What the algorithm is:
# QuickSelect is a randomized selection algorithm that finds the k-th
# smallest element in expected linear time. It's related to Quick Sort
# but only recurses into one partition instead of both.
#
# How it works conceptually:
# Pick a random pivot element. Partition the array into elements less
# than the pivot, equal to the pivot, and greater than the pivot. If
# k falls in the "less than" partition, recurse there. If k falls in
# the "equal" partition, return the pivot. Otherwise, recurse into the
# "greater than" partition with an adjusted k.
#
# Time and space complexity:
# Time: O(n) expected, O(n^2) worst case.
# Space: O(1) — in-place partitioning.
#
# Where this appears in bioinformatics:
# Selection algorithms find medians, quartiles, and other order statistics
# in genomic data — for example, finding the median expression level of
# a gene across samples, or the 95th percentile of read quality scores.


import random
import sys
import os


def quickselect(arr, k):
    """Find the k-th smallest element (1-indexed) using QuickSelect.

    Args:
        arr: A list of integers.
        k: The rank (1-indexed) of the element to find.

    Returns:
        The k-th smallest element.
    """
    if len(arr) == 1:
        return arr[0]

    pivot = random.choice(arr)

    lows = [x for x in arr if x < pivot]
    highs = [x for x in arr if x > pivot]
    pivots = [x for x in arr if x == pivot]

    if k <= len(lows):
        return quickselect(lows, k)
    elif k <= len(lows) + len(pivots):
        return pivots[0]
    else:
        return quickselect(highs, k - len(lows) - len(pivots))


def main():
    input_file = 'med_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            content = f.read()
    else:
        content = sys.stdin.read()

    lines = content.strip().splitlines()
    n = int(lines[0].strip())
    arr = list(map(int, lines[1].strip().split()))
    k = int(lines[2].strip())

    result = quickselect(arr, k)
    print(result)


if __name__ == '__main__':
    main()
