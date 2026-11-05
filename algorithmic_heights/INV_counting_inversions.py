# INV (Counting Inversions)
# https://rosalind.info/problems/inv/
#
# What the algorithm is:
# An inversion is a pair (i, j) with i < j and A[i] > A[j]. Counting
# inversions measures how far an array is from being sorted. The
# efficient approach uses merge sort: during the merge step, whenever
# we take an element from the right half before all elements of the
# left half, that element forms an inversion with every remaining
# element in the left half.
#
# How it works conceptually:
# Split the array into halves, recursively count inversions in each
# half, then count inversions that cross the midpoint during the merge.
# When merging, if we pick an element from the right half while elements
# remain in the left half, each remaining left element forms an inversion
# with this right element (since all left elements are greater than this
# right element and appear before it).
#
# Time and space complexity:
# Time: O(n log n) — same as merge sort, with O(1) extra work per merge
# step.
# Space: O(n) — temporary arrays for merging.
#
# Where this appears in bioinformatics:
# Inversion counting is used in genome rearrangement studies — for
# example, measuring the number of inversions between two gene order
# permutations, which relates to reversal distance in evolutionary
# analysis.


def merge_count(arr, temp, left, mid, right):
    """Merge two sorted halves and count cross-inversions.

    Args:
        arr: The array being sorted.
        temp: A temporary array for merging.
        left: Start index of left half.
        mid: End index of left half (exclusive).
        right: End index of right half (exclusive).

    Returns:
        The number of cross-inversions.
    """
    i = left
    j = mid
    k = left
    inversions = 0

    while i < mid and j < right:
        if arr[i] <= arr[j]:
            temp[k] = arr[i]
            i += 1
        else:
            temp[k] = arr[j]
            j += 1
            # All remaining elements in left half form inversions
            # with arr[j] (which was originally at position j)
            inversions += mid - i
        k += 1

    # Copy remaining elements
    while i < mid:
        temp[k] = arr[i]
        i += 1
        k += 1
    while j < right:
        temp[k] = arr[j]
        j += 1
        k += 1

    # Copy back to original array
    for i in range(left, right):
        arr[i] = temp[i]

    return inversions


def merge_sort_count(arr, temp, left, right):
    """Sort arr[left:right] and count inversions using merge sort.

    Args:
        arr: The array to sort.
        temp: A temporary array.
        left: Start index (inclusive).
        right: End index (exclusive).

    Returns:
        The number of inversions in arr[left:right].
    """
    if right - left <= 1:
        return 0

    mid = (left + right) // 2
    inversions = 0

    inversions += merge_sort_count(arr, temp, left, mid)
    inversions += merge_sort_count(arr, temp, mid, right)
    inversions += merge_count(arr, temp, left, mid, right)

    return inversions


def count_inversions(arr):
    """Count the number of inversions in arr.

    Args:
        arr: A list of integers.

    Returns:
        The number of pairs (i, j) with i < j and arr[i] > arr[j].
    """
    temp = [0] * len(arr)
    # Work on a copy to avoid modifying the original
    arr_copy = arr[:]
    return merge_sort_count(arr_copy, temp, 0, len(arr_copy))


def main():
    import sys
    import os

    input_file = 'inv_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            content = f.read()
    else:
        content = sys.stdin.read()

    lines = content.strip().splitlines()
    n = int(lines[0].strip())
    arr = list(map(int, lines[1].strip().split()))

    result = count_inversions(arr)
    print(result)


if __name__ == '__main__':
    main()
