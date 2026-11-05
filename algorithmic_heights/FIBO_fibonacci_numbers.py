# FIBO (Fibonacci Numbers)
# https://rosalind.info/problems/fibo/
#
# What the algorithm is:
# The Fibonacci sequence is defined by the recurrence F(n) = F(n-1) + F(n-2)
# with base cases F(0) = 0 and F(1) = 1. Each number is the sum of the two
# preceding ones. The sequence appears throughout computer science and
# mathematics — in algorithm analysis, data structures (Fibonacci heaps),
# and even in modeling biological population growth.
#
# How it works conceptually:
# The naive recursive implementation recomputes the same values over and
# over, leading to exponential time. The efficient approach uses dynamic
# programming: compute F(0), F(1), ..., F(n) iteratively, storing each
# result as you go. Each value is computed exactly once, giving linear
# time O(n) and constant space O(1) if you only keep the last two values.
#
# Time and space complexity:
# Time: O(n) — a single loop from 2 to n, each iteration doing constant work.
# Space: O(1) — only the two most recent Fibonacci numbers need to be stored.
#
# Where this appears in bioinformatics:
# Fibonacci numbers model idealized population growth (the original rabbit
# problem). They also appear in the analysis of Euclid's algorithm for
# computing GCD, where consecutive Fibonacci numbers are the worst-case
# input. In bioinformatics, Fibonacci-like recurrences appear in population
# genetics and in modeling the growth of branching structures.


def fibonacci(n):
    """Compute the nth Fibonacci number using iterative dynamic programming.

    F(0) = 0, F(1) = 1, F(n) = F(n-1) + F(n-2) for n > 1.

    Args:
        n: A non-negative integer (n <= 25 for this problem).

    Returns:
        The nth Fibonacci number.
    """
    if n == 0:
        return 0
    if n == 1:
        return 1

    prev, curr = 0, 1
    for _ in range(2, n + 1):
        prev, curr = curr, prev + curr

    return curr


def main():
    import sys
    import os

    input_file = 'fibo_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            n = int(f.read().strip())
    else:
        n = int(sys.stdin.read().strip())

    print(fibonacci(n))


if __name__ == '__main__':
    main()
