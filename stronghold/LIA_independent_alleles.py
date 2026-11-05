# LIA (Independent Alleles)
# https://rosalind.info/problems/lia/
#
# What independent alleles means:
# Mendel's second law, the law of independent assortment, says that alleles
# for different traits are inherited independently of each other. When you
# cross two organisms that are heterozygous for two factors (Aa Bb × Aa Bb),
# the inheritance of the A/a factor doesn't affect the inheritance of the
# B/b factor. This lets us multiply probabilities: P(Aa Bb) = P(Aa) × P(Bb).
#
# Why this matters:
# Independent assortment is why traits can be combined in new ways in
# offspring. It's the basis for predicting complex genetic outcomes and for
# understanding how genetic variation is generated in populations. Without
# independent assortment, all genes on the same chromosome would always be
# inherited together, and evolution would have much less raw material to
# work with.
#
# What this solution does:
# Each child of an Aa Bb × Aa Bb cross has probability 1/4 of being Aa Bb
# (1/2 for Aa × 1/2 for Bb). Generation k has 2^k organisms (each organism
# has 2 children, starting from Tom in generation 0). The number of Aa Bb
# organisms in generation k follows a binomial distribution with 2^k trials
# and success probability 1/4. We compute P(X ≥ N) = 1 - P(X < N).


def independent_alleles(k, n):
    """Return the probability that at least n Aa Bb organisms are in
    generation k of Tom's family tree.

    Each organism in generation k has probability 1/4 of being Aa Bb.
    Generation k has 2^k organisms. X ~ Binomial(2^k, 1/4).

    Args:
        k: Generation number (k >= 0).
        n: Minimum number of Aa Bb organisms required.

    Returns:
        Probability that at least n organisms in generation k are Aa Bb.
    """
    import math

    total_organisms = 2 ** k
    p = 0.25  # P(Aa Bb) for each child

    # P(X >= n) = 1 - P(X < n) = 1 - sum_{i=0}^{n-1} C(total, i) * p^i * (1-p)^(total-i)
    prob_less_than_n = 0.0
    for i in range(n):
        prob_less_than_n += (
            math.comb(total_organisms, i)
            * (p ** i)
            * ((1 - p) ** (total_organisms - i))
        )

    return 1.0 - prob_less_than_n


def main():
    import sys
    import os

    input_file = 'lia_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            k, n = map(int, f.read().split())
    else:
        k, n = map(int, sys.stdin.read().split())

    result = independent_alleles(k, n)
    print(f"{result:.3f}")


if __name__ == '__main__':
    main()
