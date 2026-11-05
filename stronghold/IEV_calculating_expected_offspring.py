# IEV (Calculating Expected Offspring)
# https://rosalind.info/problems/iev/
#
# What this problem is about:
# In genetics, when two organisms mate, their offspring inherit one allele
# from each parent. The combination of alleles determines the offspring's
# genotype (AA, Aa, or aa), and the genotype determines the phenotype
# (whether the dominant trait shows). A dominant allele (capital letter,
# like A) expresses its trait even with just one copy. A recessive allele
# (lowercase, like a) only shows when two copies are present (aa).
#
# Why expected value matters:
# Genetic crosses are random — you can't predict the exact genotype of any
# single offspring. But you can predict the average outcome over many
# offspring. The expected value tells you the long-term average number of
# offspring with the dominant phenotype. This is essential in genetic
# counseling: if both parents are carriers of a recessive condition, what's
# the expected number of affected children?
#
# The six genotype pairings and their dominant offspring probabilities:
# 1. AA × AA: All offspring are AA. Dominant probability = 1.
# 2. AA × Aa: Offspring are 1/2 AA, 1/2 Aa. All dominant. Probability = 1.
# 3. AA × aa: All offspring are Aa. All dominant. Probability = 1.
# 4. Aa × Aa: Offspring are 1/4 AA, 1/2 Aa, 1/4 aa. Dominant = 3/4.
# 5. Aa × aa: Offspring are 1/2 Aa, 1/2 aa. Dominant = 1/2.
# 6. aa × aa: All offspring are aa. Dominant probability = 0.
#
# Each couple has exactly 2 offspring, so multiply each probability by 2
# and by the number of couples of that type.


def expected_offspring(counts):
    """Calculate expected number of offspring with dominant phenotype.

    Args:
        counts: A list of 6 integers representing couples with genotype
                pairings in order: AA-AA, AA-Aa, AA-aa, Aa-Aa, Aa-aa, aa-aa.

    Returns:
        The expected number of offspring displaying the dominant phenotype.
    """
    # Expected dominant offspring per couple for each genotype pairing
    dominant_prob = [1.0, 1.0, 1.0, 0.75, 0.5, 0.0]
    offspring_per_couple = 2

    total = 0.0
    for i in range(6):
        total += counts[i] * dominant_prob[i] * offspring_per_couple

    return total


def main():
    import sys
    import os

    input_file = 'iev_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            counts = list(map(int, f.read().split()))
    else:
        counts = list(map(int, sys.stdin.read().split()))

    result = expected_offspring(counts)
    print(result)


if __name__ == '__main__':
    main()
