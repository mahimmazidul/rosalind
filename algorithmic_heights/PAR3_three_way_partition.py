# PAR3 (Three Way Partition)
# https://rosalind.info/problems/par3/
#
# What the algorithm is:
# Three-way partition (also known as Dutch National Flag) rearranges an
# array around a pivot element so that all elements less than the pivot
# come first, then all elements equal to the pivot, then all elements
# greater than the pivot. This is the partitioning scheme used in
# 3-way Quick Sort, which is efficient when there are many duplicate
# keys.
#
# How it works conceptually:
# Use three pointers: low (boundary of elements < pivot), mid (current
# element being examined), and high (boundary of elements > pivot).
# - If arr[mid] < pivot: swap arr[low] and arr[mid], advance both.
# - If arr[mid] == pivot: advance mid.
# - If arr[mid] > pivot: swap arr[mid] and arr[high], decrement high.
#
# Time and space complexity:
# Time: O(n) — single pass through the array.
# Space: O(1) — in-place, only a few integer variables.
#
# Where this appears in bioinformatics:
# 3-way partitioning is used in Quick Sort variants for sorting genomic
# data with many duplicate values (e.g., sorted by chromosome, where
# many positions share the same chromosome number).


def three_way_partition(arr):
    """Partition arr around the first element using Dutch National Flag.

    Elements < pivot come first, then elements == pivot, then elements
    > pivot.

    Args:
        arr: A list of integers.

    Returns:
        The permuted list with elements grouped around the pivot.
    """
    if len(arr) <= 1:
        return arr

    pivot = arr[0]
    low = 0
    mid = 0
    high = len(arr) - 1

    while mid <= high:
        if arr[mid] < pivot:
            arr[low], arr[mid] = arr[mid], arr[low]
            low += 1
            mid += 1
        elif arr[mid] == pivot:
            mid += 1
        else:
            arr[mid], arr[high] = arr[high], arr[mid]
            high -= 1

    return arr


def main():
    import sys
    import os

    input_file = 'par3_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            content = f.read()
    else:
        content = sys.stdin.read()

    lines = content.strip().splitlines()
    n = int(lines[0].strip())
    arr = list(map(int, lines[1].strip().split()))

    result = three_way_partition(arr)
    print(' '.join(str(x) for x in result))


if __name__ == '__main__':
    main()
