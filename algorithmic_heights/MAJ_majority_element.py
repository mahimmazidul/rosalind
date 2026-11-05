# MAJ (Majority Element)
# https://rosalind.info/problems/maj/
#
# What the algorithm is:
# The majority element of an array is an element that appears more than
# half the time. The Boyer-Moore majority vote algorithm finds a candidate
# in a single pass using constant space, then verifies it with a second
# pass.
#
# How it works conceptually:
# Maintain a candidate element and a counter. Start with the first element
# as candidate and counter = 1. For each subsequent element: if it equals
# the candidate, increment the counter; otherwise, decrement. If the
# counter reaches zero, switch the candidate to the current element and
# reset the counter to 1. The intuition: pairs of different elements
# cancel each other out, so if a majority element exists, it will survive
# as the final candidate.
#
# Time and space complexity:
# Time: O(n) — two linear passes (one to find candidate, one to verify).
# Space: O(1) — only the candidate and counter variables.
#
# Where this appears in bioinformatics pipelines:
# Majority voting appears in consensus building — for example, determining
# the most common base at a position across multiple aligned sequencing
# reads, or consensus calling in variant detection where multiple reads
# support the same allele.


def majority_element(arr):
    """Find the majority element using Boyer-Moore majority vote algorithm.

    Returns the element that appears more than n/2 times, or -1 if none
    exists.

    Args:
        arr: A list of integers.

    Returns:
        The majority element, or -1 if no majority element exists.
    """
    # First pass: find candidate
    candidate = None
    count = 0

    for num in arr:
        if count == 0:
            candidate = num
            count = 1
        elif num == candidate:
            count += 1
        else:
            count -= 1

    # Second pass: verify the candidate
    actual_count = sum(1 for num in arr if num == candidate)
    n = len(arr)

    if actual_count > n // 2:
        return candidate
    return -1


def main():
    import sys
    import os

    input_file = 'maj_input.txt'
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
        results.append(str(majority_element(arr)))

    print(' '.join(results))


if __name__ == '__main__':
    main()
