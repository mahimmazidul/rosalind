# BINS (Binary Search)
# https://rosalind.info/problems/bins/
#
# What the algorithm is:
# Binary search is a divide-and-conquer algorithm for finding a target value
# in a sorted array. It repeatedly compares the target with the middle element
# of the current search range and discards half of the remaining elements each
# time. If the target equals the middle element, the search succeeds. If the
# target is smaller, search the left half; if larger, search the right half.
#
# How it works conceptually:
# Think of looking up a name in a phone book: you open to the middle, decide
# whether your target is before or after that page, and repeat with the
# relevant half. Each step halves the search space, so the number of steps
# grows logarithmically with the array size.
#
# Time and space complexity:
# Time: O(log n) per query — each comparison halves the search range.
# Space: O(1) — only a few integer variables for the search bounds.
# With m queries, total time is O(m log n).
#
# Where this appears in bioinformatics pipelines:
# Binary search is used whenever you need to find a specific entry in a
# sorted index — for example, looking up a gene by name in a sorted
# annotation file, finding a position in a sorted list of variants, or
# querying a suffix array for pattern matches. It's a fundamental building
# block of many bioinformatics tools that work with large, sorted datasets.


def binary_search(arr, target):
    """Search for target in sorted array arr using binary search.

    Returns the 1-based index of target if found, or -1 if not found.

    Args:
        arr: A sorted list of integers.
        target: The integer value to search for.

    Returns:
        The 1-based index of target in arr, or -1 if not present.
    """
    left = 0
    right = len(arr) - 1

    while left <= right:
        mid = (left + right) // 2
        if arr[mid] == target:
            return mid + 1  # 1-based index
        elif arr[mid] < target:
            left = mid + 1
        else:
            right = mid - 1

    return -1


def main():
    import sys
    import os

    input_file = 'bins_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            content = f.read()
    else:
        content = sys.stdin.read()

    lines = content.strip().splitlines()
    n = int(lines[0].strip())
    m = int(lines[1].strip())
    arr = list(map(int, lines[2].strip().split()))
    queries = list(map(int, lines[3].strip().split()))

    results = [str(binary_search(arr, q)) for q in queries]
    print(' '.join(results))


if __name__ == '__main__':
    main()
