# 2SUM (Two Sum)
# https://rosalind.info/problems/2sum/
#
# What the algorithm is:
# Given an array, find two indices p < q such that A[p] = -A[q], i.e.,
# the two values sum to zero. This is a classic hash table lookup problem:
# for each element, check if its negation has already been seen.
#
# How it works conceptually:
# Iterate through the array. For each element x at position i (1-based),
# check if -x was already seen at some earlier position j. If so, return
# (j, i) since j < i. Otherwise, record x at position i and continue.
# Using a hash map gives O(1) average lookup per element.
#
# Time and space complexity:
# Time: O(n) per array — single pass with O(1) hash lookups.
# Space: O(n) — hash map stores up to n elements.
# For k arrays, total time is O(k * n).
#
# Where this appears in bioinformatics:
# Two-sum style lookups appear when searching for complementary values
# — for example, finding pairs of mutations that cancel each other's
# effect, or matching mass spectrometry peaks that sum to a known
# precursor mass.


def two_sum(arr):
    """Find two indices p < q such that arr[p] = -arr[q].

    Uses a hash map to store the first position of each value seen.

    Args:
        arr: A list of integers.

    Returns:
        A tuple (p, q) of 1-based indices, or -1 if no such pair exists.
    """
    seen = {}  # value -> 1-based index

    for i, x in enumerate(arr):
        target = -x
        if target in seen:
            return (seen[target], i + 1)
        if x not in seen:
            seen[x] = i + 1

    return -1


def main():
    import sys
    import os

    input_file = '2sum_input.txt'
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
        result = two_sum(arr)
        if result == -1:
            results.append('-1')
        else:
            results.append(f'{result[0]} {result[1]}')

    print(' '.join(results))


if __name__ == '__main__':
    main()
