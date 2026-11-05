# HEA (Building a Heap)
# https://rosalind.info/problems/hea/
#
# What the algorithm is:
# A binary heap is a complete binary tree stored in an array, where each
# parent node is greater than or equal to its children (max heap). Building
# a heap from an arbitrary array means rearranging elements so the max heap
# property holds: for every index i > 0, A[parent(i)] >= A[i].
#
# How it works conceptually:
# The problem specifies a particular approach: go from the end of the array
# to the beginning, and for each element, "bubble it up" by repeatedly
# swapping it with its parent if it's larger, until it reaches a position
# where its parent is larger or it becomes the root. This is essentially
# inserting each element into a heap one at a time, processing in reverse
# order.
#
# The parent of node j (1-indexed) is ceil(j/2), and children are 2j and
# 2j+1. This defines the tree structure used for bubbling.
#
# Time and space complexity:
# Time: O(n log n) — each of n elements bubbles up at most O(log n) levels.
# (A tighter analysis shows O(n), but the simple bound is O(n log n).)
# Space: O(1) — in-place, only a few variables for swapping.
#
# Where this appears in bioinformatics pipelines:
# Heaps implement priority queues, used in Dijkstra's shortest path (for
# sequence alignment graphs), Prim's algorithm (for minimum spanning trees
# in phylogenetic networks), and event-driven simulations in population
# genetics.


def parent_index(j):
    """Return the 0-indexed parent index of 0-indexed position j.

    Using the problem's convention: 1-indexed parent of j is ceil(j/2).
    Converting to 0-indexed: parent = ceil((j+1)/2) - 1 = (j+2)//2 - 1.
    """
    return (j + 2) // 2 - 1


def build_max_heap(arr):
    """Build a max heap from arr using bubble-up from end to beginning.

    Args:
        arr: A list of integers.

    Returns:
        The same list rearranged to satisfy the max heap property.
    """
    n = len(arr)
    # Go from end to beginning (0-indexed: n-1 down to 1)
    for i in range(n - 1, 0, -1):
        # Bubble up element at position i
        j = i
        while j > 0:
            p = parent_index(j)
            if arr[p] < arr[j]:
                arr[p], arr[j] = arr[j], arr[p]
                j = p
            else:
                break
    return arr


def main():
    import sys
    import os

    input_file = 'hea_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            content = f.read()
    else:
        content = sys.stdin.read()

    lines = content.strip().splitlines()
    n = int(lines[0].strip())
    arr = list(map(int, lines[1].strip().split()))

    result = build_max_heap(arr)
    print(' '.join(str(x) for x in result))


if __name__ == '__main__':
    main()
