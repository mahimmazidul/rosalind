# MS (Merge Sort)
# https://rosalind.info/problems/ms/
#
# What the algorithm is:
# Merge sort is a divide-and-conquer sorting algorithm. It splits the array
# into two halves, recursively sorts each half, then merges the two sorted
# halves into a single sorted array.
#
# How it works conceptually:
# 1. If the array has 0 or 1 elements, it's already sorted.
# 2. Split the array into two halves.
# 3. Recursively sort each half.
# 4. Merge the two sorted halves using the two-pointer technique.
#
# Time and space complexity:
# Time: O(n log n) — the recurrence T(n) = 2T(n/2) + O(n) solves to O(n log n).
# Space: O(n) — temporary arrays for merging.
#
# Where this appears in bioinformatics pipelines:
# Merge sort is used when stable sorting of large genomic datasets is
# needed — for example, sorting millions of variants by genomic position,
# or sorting sequencing reads by alignment coordinate. Its O(n log n)
# guarantee makes it predictable for large inputs.


def merge(left, right):
    """Merge two sorted lists into a single sorted list."""
    result = []
    i, j = 0, 0
    while i < len(left) and j < len(right):
        if left[i] <= right[j]:
            result.append(left[i])
            i += 1
        else:
            result.append(right[j])
            j += 1
    result.extend(left[i:])
    result.extend(right[j:])
    return result


def merge_sort(arr):
    """Sort arr using merge sort.

    Args:
        arr: A list of integers.

    Returns:
        A new sorted list.
    """
    if len(arr) <= 1:
        return arr

    mid = len(arr) // 2
    left = merge_sort(arr[:mid])
    right = merge_sort(arr[mid:])

    return merge(left, right)


def main():
    import sys
    import os

    input_file = 'ms_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            content = f.read()
    else:
        content = sys.stdin.read()

    lines = content.strip().splitlines()
    n = int(lines[0].strip())
    arr = list(map(int, lines[1].strip().split()))

    sorted_arr = merge_sort(arr)
    print(' '.join(str(x) for x in sorted_arr))


if __name__ == '__main__':
    main()
