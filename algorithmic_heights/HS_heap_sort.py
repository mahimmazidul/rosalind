# HS (Heap Sort)
# https://rosalind.info/problems/hs/
#
# What the algorithm is:
# Heap sort sorts an array by first transforming it into a max heap, then
# repeatedly extracting the maximum element (the root) and placing it at
# the end of the array.
#
# How it works conceptually:
# Phase 1: Build a max heap from the array — the largest element is at
# the root. Phase 2: Swap the root with the last element, reduce the heap
# size by 1, and sift-down the new root to restore the heap property.
# Repeat until the heap is empty. The array ends up sorted in ascending
# order.
#
# Time and space complexity:
# Time: O(n log n) — building the heap is O(n), each of n-1 sift-down
# operations is O(log n).
# Space: O(1) — in-place sorting.
#
# Graph/heap representation:
# The array represents a complete binary tree where node i (0-indexed)
# has parent (i-1)//2, left child 2i+1, and right child 2i+2.
#
# Where this appears in bioinformatics:
# Heap sort provides guaranteed O(n log n) sorting for large genomic
# datasets where worst-case performance matters more than cache
# efficiency (unlike quick sort which can degrade to O(n^2)).


def sift_down(arr, start, end):
    """Sift down the element at start to restore max heap property.

    The heap occupies indices 0 through end-1.

    Args:
        arr: The array representing the heap.
        start: The index of the element to sift down.
        end: One past the last valid heap index.
    """
    root = start
    while True:
        left = 2 * root + 1
        right = 2 * root + 2
        largest = root

        if left < end and arr[left] > arr[largest]:
            largest = left
        if right < end and arr[right] > arr[largest]:
            largest = right

        if largest == root:
            break

        arr[root], arr[largest] = arr[largest], arr[root]
        root = largest


def build_max_heap(arr):
    """Transform arr into a max heap in-place using sift-down."""
    n = len(arr)
    # Start from the last non-leaf node and sift down each
    for i in range(n // 2 - 1, -1, -1):
        sift_down(arr, i, n)


def heap_sort(arr):
    """Sort arr using heap sort (in-place).

    Args:
        arr: A list of integers.

    Returns:
        The same list, sorted in ascending order.
    """
    n = len(arr)

    # Phase 1: Build max heap
    build_max_heap(arr)

    # Phase 2: Repeatedly extract max and place at end
    for end in range(n - 1, 0, -1):
        # Swap root (max) with last element
        arr[0], arr[end] = arr[end], arr[0]
        # Sift down the new root in the reduced heap
        sift_down(arr, 0, end)

    return arr


def main():
    import sys
    import os

    input_file = 'hs_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            content = f.read()
    else:
        content = sys.stdin.read()

    lines = content.strip().splitlines()
    n = int(lines[0].strip())
    arr = list(map(int, lines[1].strip().split()))

    heap_sort(arr)
    print(' '.join(str(x) for x in arr))


if __name__ == '__main__':
    main()
