# IPRB (Mendel's First Law)
# https://rosalind.info/problems/iprb/
#
# Mendel's first law — the law of segregation:
# Every organism carries two copies (alleles) of each gene, one inherited
# from each parent. When an organism produces gametes (sperm or egg cells),
# the two alleles segregate — they separate — so each gamete gets only one
# of the two. During fertilization, offspring receive one allele from each
# parent, restoring the pair. This is why traits can skip generations and
# reappear: a recessive allele can be carried silently and passed on without
# being expressed.
#
# Dominant and recessive alleles:
# A dominant allele (written as a capital letter, like Y) expresses its trait
# even when only one copy is present. A recessive allele (lowercase, like y)
# only expresses when two copies are present. An organism with genotype YY
# is homozygous dominant, Yy is heterozygous, and yy is homozygous recessive.
# Heterozygotes (Yy) look dominant but carry the recessive allele and can
# pass it to their offspring. This is why we care about dominant probability
# in genetic counseling: two seemingly unaffected parents can have a child
# with a recessive condition if both are carriers.
#
# Punnett square logic, generalised to populations:
# A Punnett square shows all possible outcomes of one mating pair. This
# problem generalises that idea: instead of one pair, we sample two
# individuals from a whole population and ask what's the probability their
# child shows the dominant trait. We enumerate every possible pairing of
# genotypes and weight each by how likely that pairing is to occur.
#
# Why we sample without replacement:
# We're choosing two mating individuals from a finite population. Once the
# first individual is chosen, they can't be chosen again as the second
# parent (an organism doesn't mate with itself in this model). So the total
# number of ordered pairs is N*(N-1), where N = k + m + n, and the
# probabilities shift after the first draw.
#
# Real applications:
# - Predicting inheritance of genetic diseases: if both parents are carriers
#   of a recessive condition (like cystic fibrosis), each child has a 25%
#   chance of being affected.
# - Understanding why harmful recessive diseases persist: carriers are
#   unaffected, so natural selection doesn't remove the allele from the
#   population; it hides in heterozygotes and gets passed on silently.
#
# This is Rosalind problem #6. Given k homozygous dominant, m heterozygous,
# and n homozygous recessive individuals, compute the probability that two
# randomly selected mating individuals produce an offspring with at least
# one dominant allele.


def dominant_probability(k, m, n):
    """Return the probability of a dominant-allele offspring.

    We compute P(dominant offspring) = 1 - P(recessive offspring).
    An offspring is recessive (yy) only if both parents pass on a y allele.
    Since we select two distinct individuals without replacement, we enumerate
    every ordered pair of genotypes and weight by the probability of that
    pairing occurring, then by the probability both parents contribute y.

    Args:
        k: Number of homozygous dominant individuals (YY).
        m: Number of heterozygous individuals (Yy).
        n: Number of homozygous recessive individuals (yy).

    Returns:
        Probability that a random mating produces an offspring with at
        least one dominant allele.
    """
    total = k + m + n

    # Total number of ordered pairs of distinct individuals
    total_pairs = total * (total - 1)

    # Count ordered pairs where both parents can pass on a y allele,
    # weighted by the probability both actually do:

    # Yy x Yy: m*(m-1) ordered pairs, each passes y with probability 1/2 * 1/2 = 1/4
    yy_from_heterozygous_pairs = m * (m - 1) * 0.25

    # Yy x yy and yy x Yy: m*n + n*m ordered pairs, each passes y with
    # probability 1/2 * 1 = 1/2
    yy_from_mixed_pairs = m * n * 0.5 + n * m * 0.5

    # yy x yy: n*(n-1) ordered pairs, each passes y with probability 1 * 1 = 1
    yy_from_recessive_pairs = n * (n - 1) * 1.0

    # YY cannot pass on y, so pairs involving YY never produce yy offspring.

    # P(both parents pass y) = weighted count / total ordered pairs
    p_recessive_offspring = (
        yy_from_heterozygous_pairs + yy_from_mixed_pairs + yy_from_recessive_pairs
    ) / total_pairs

    # P(at least one dominant allele) = 1 - P(recessive offspring)
    return 1.0 - p_recessive_offspring


def main():
    import sys
    import os

    input_file = 'iprb_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            k, m, n = map(int, f.read().split())
    else:
        k, m, n = map(int, sys.stdin.read().split())

    result = dominant_probability(k, m, n)
    print(f"{result:.5f}")


if __name__ == '__main__':
    main()
