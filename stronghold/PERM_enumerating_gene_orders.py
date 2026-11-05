# PERM (Enumerating Gene Orders)
# https://rosalind.info/problems/perm/
#
# What a permutation is:
# A permutation is an ordering of a set of distinct items. In this problem,
# the items are the integers 1 through n, and a permutation of length n is
# any arrangement of those n numbers in a sequence. For example, (3, 1, 2)
# is one permutation of length 3.
#
# Why permutations matter in genome rearrangements:
# When comparing two genomes, researchers identify synteny blocks — similar
# intervals of DNA — and label each with a number. A genome then becomes a
# permutation of these numbered blocks. Genome rearrangements (like
# reversals, transpositions, and translocations) change the order of these
# blocks, turning one permutation into another. Understanding the space of
# all possible permutations helps us measure how far apart two genomes are
# in terms of rearrangement events.
#
# Why enumerate all permutations:
# For small n (up to 7), we can list every possible gene order. This is
# useful for understanding the rearrangement space, testing algorithms, and
# computing distances between genomes. The total number of permutations is
# n! (n factorial), which grows quickly: 7! = 5040, still manageable.


import itertools


def enumerate_permutations(n):
    """Generate all permutations of the integers 1 through n.

    Args:
        n: A positive integer (n <= 7).

    Returns:
        A tuple of (count, permutations) where count is n! and
        permutations is a list of tuples, each a permutation of 1..n.
    """
    perms = list(itertools.permutations(range(1, n + 1)))
    return len(perms), perms


def main():
    import sys
    import os

    input_file = 'perm_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            n = int(f.read().strip())
    else:
        n = int(sys.stdin.read().strip())

    count, perms = enumerate_permutations(n)

    print(count)
    for perm in perms:
        print(' '.join(str(x) for x in perm))


if __name__ == '__main__':
    main()
