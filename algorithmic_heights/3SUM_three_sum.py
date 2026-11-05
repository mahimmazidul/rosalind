# 3SUM (Three Sum)
# https://rosalind.info/problems/3sum/
#
# What the algorithm is:
# 3SUM asks: given an array, find three distinct elements that sum to zero.
# It's a classic problem in computational geometry and algorithm design.
#
# How it works conceptually:
# Sort the array while tracking original indices. For each element at
# position i, use two pointers (left from i+1, right from end) to find
# two more elements that sum to -arr[i]. Move the pointers based on
# whether the current sum is too small or too large.
#
# Time and space complexity:
# Time: O(n^2) — sorting is O(n log n), then for each of n elements,
# the two-pointer scan is O(n).
# Space: O(n) — for the indexed array.
#
# Where this appears in bioinformatics:
# 3SUM-like problems appear in mass spectrometry — finding three peptide
# fragments whose masses sum to a target, or in population genetics where
# three allele frequencies might sum to a known total.


def three_sum(arr):
    """Find three indices p < q < r such that arr[p] + arr[q] + arr[r] = 0.

    Args:
        arr: A list of integers.

    Returns:
        A tuple of three 1-based indices (p, q, r) with p < q < r,
        or -1 if no such triple exists.
    """
    n = len(arr)
    # Sort by value, keeping original 1-based indices
    indexed = sorted([(arr[i], i + 1) for i in range(n)], key=lambda x: x[0])

    for i in range(n - 2):
        left = i + 1
        right = n - 1
        while left < right:
            s = indexed[i][0] + indexed[left][0] + indexed[right][0]
            if s == 0:
                result = sorted([indexed[i][1], indexed[left][1], indexed[right][1]])
                return tuple(result)
            elif s < 0:
                left += 1
            else:
                right -= 1

    return -1


def main():
    import sys
    import os

    input_file = '3sum_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            content = f.read()
    else:
        content = sys.stdin.read()

    lines = content.strip().splitlines()
    k, n = map(int, lines[0].split())

    results = []
    for i in range(1, k + 1):
        arr = list(map(int, lines[i].strip().split()))
        result = three_sum(arr)
        if result == -1:
            results.append('-1')
        else:
            results.append(f'{result[0]} {result[1]} {result[2]}')

    print(' '.join(results))


if __name__ == '__main__':
    main()
