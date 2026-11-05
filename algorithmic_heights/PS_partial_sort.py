# PS (Partial Sort)
# https://rosalind.info/problems/ps/
#
# What the algorithm is:
# Partial sort returns the k smallest elements of an array in sorted
# order. Since k ≤ 1000 (much smaller than n ≤ 10^5), we can simply
# sort the entire array and take the first k elements.
#
# How it works conceptually:
# Sort the array in ascending order, then return the first k elements.
# These are guaranteed to be the k smallest elements in sorted order.
#
# Time and space complexity:
# Time: O(n log n) — sorting the full array.
# Space: O(n) — the sorted array.
#
# A more efficient approach for small k would be O(n log k) using a
# max-heap of size k, but O(n log n) is fast enough for n ≤ 10^5.
#
# Where this appears in bioinformatics:
# Partial sort is used when you only need the top-k results — for
# example, the top-k most differentially expressed genes, or the
# k highest-scoring sequence alignments.


def main():
    import sys
    import os

    input_file = 'ps_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            content = f.read()
    else:
        content = sys.stdin.read()

    lines = content.strip().splitlines()
    n = int(lines[0].strip())
    arr = list(map(int, lines[1].strip().split()))
    k = int(lines[2].strip())

    arr.sort()
    result = arr[:k]
    print(' '.join(str(x) for x in result))


if __name__ == '__main__':
    main()
