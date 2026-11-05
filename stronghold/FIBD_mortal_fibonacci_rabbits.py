# FIBD (Mortal Fibonacci Rabbits)
# https://rosalind.info/problems/fibd/
#
# Fibonacci's rabbit model and its origin:
# In 1202, Leonardo of Pisa (Fibonacci) introduced what we now call the
# Fibonacci sequence through a thought experiment about rabbit breeding.
# The model assumes a newborn pair of rabbits matures in one month and then
# produces a new pair every month thereafter. Rabbits never die in the
# original model, so the population grows as F(n) = F(n-1) + F(n-2).
# This was never meant to be a realistic population model — it was a
# mathematical puzzle about idealized breeding.
#
# Why we need mortal rabbits:
# Real rabbits don't live forever. In this version, every rabbit dies after
# a fixed number of months m. This makes the model more biologically
# plausible and turns it into a proper population dynamics problem. Instead
# of unbounded exponential growth, the population reaches a balance between
# births and deaths. This is the kind of model population geneticists and
# ecologists use when studying how fast a population can grow given finite
# lifespans and reproductive rates.
#
# What this solution does:
# It tracks rabbit pairs by age (newborn, 1 month old, ..., m-1 months old).
# Each month, mature rabbits (age 1 or older) produce one newborn pair.
# All rabbits age by one month, and any rabbit that reaches age m dies.
# After n months, the total number of living pairs is returned.


def mortal_fibonacci(n, m):
    """Return the number of rabbit pairs alive after n months with lifespan m.

    Uses dynamic programming by tracking the population by age. Each month,
    mature rabbits (age >= 1) produce newborns, all rabbits age by one, and
    rabbits reaching age m die.

    Args:
        n: Number of months to simulate.
        m: Lifespan of a rabbit in months.

    Returns:
        The total number of rabbit pairs alive after month n.
    """
    # pop[age] = number of rabbit pairs at that age
    # age 0 = newborn, age m-1 = oldest (dies after this month)
    pop = [0] * m
    pop[0] = 1  # start with one newborn pair in month 1

    for _ in range(2, n + 1):
        newborns = sum(pop[1:])  # mature rabbits (age 1+) reproduce
        # Age all rabbits by one month; rabbits at age m-1 die
        for age in range(m - 1, 0, -1):
            pop[age] = pop[age - 1]
        pop[0] = newborns

    return sum(pop)


def main():
    import sys
    import os

    input_file = 'fibd_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            n, m = map(int, f.read().split())
    else:
        n, m = map(int, sys.stdin.read().split())

    print(mortal_fibonacci(n, m))


if __name__ == '__main__':
    main()
