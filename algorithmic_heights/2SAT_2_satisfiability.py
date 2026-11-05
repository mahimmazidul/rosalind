# 2SAT (2-Satisfiability)
# https://rosalind.info/problems/2sat/
#
# What the algorithm is:
# 2SAT solved via reduction to strongly connected components. Each clause
# (a ∨ b) implies two implications: ¬a → b and ¬b → a. Build an
# implication graph with 2n nodes (one per literal, positive and negative).
# If any SCC contains both x and ¬x, the formula is unsatisfiable.
# Otherwise, a satisfying assignment exists.
#
# How it works conceptually:
# - Map each literal to a node: x_i → 2*(i-1), ¬x_i → 2*(i-1)+1.
# - For each clause (a ∨ b), add edges ¬a → b and ¬b → a.
# - Find SCCs (Kosaraju's algorithm).
# - If comp[x_i] == comp[¬x_i] for any i, unsatisfiable.
# - Otherwise, assign x_i = TRUE if comp[x_i] > comp[¬x_i] (i.e. x_i's
#   literal is in a sink-ward SCC of the implication graph).
#
# Time and space complexity:
# Time: O(n + m) — graph construction plus two DFS passes.
# Space: O(n + m) — adjacency lists.
#
# Graph representation:
# Adjacency lists for the implication graph and its reverse.
#
# Where this appears in bioinformatics:
# 2SAT appears in haplotype inference and phylogenetics — e.g., deciding
# whether a set of pairwise compatibility constraints on genetic variants
# can be simultaneously satisfied.

import sys
import os


def solve_2sat(n, clauses):
    """Solve a 2SAT formula with n variables and m clauses.

    Args:
        n: Number of variables.
        clauses: A list of (a, b) tuples where each literal is a signed
                 integer (positive = variable, negative = negation).

    Returns:
        None if unsatisfiable, otherwise a list of signed integers
        representing the assignment (e.g. [1, -2, 3]).
    """
    N = 2 * n  # two nodes per variable

    def lit_node(lit):
        v = abs(lit)
        if lit > 0:
            return 2 * (v - 1)       # positive literal x_v
        else:
            return 2 * (v - 1) + 1   # negative literal ¬x_v

    def neg_node(node):
        return node ^ 1  # flip between x_v and ¬x_v

    # Build implication graph and its reverse
    adj = [[] for _ in range(N)]
    radj = [[] for _ in range(N)]
    for a, b in clauses:
        na, nb = lit_node(a), lit_node(b)
        # (a ∨ b) => ¬a → b and ¬b → a
        adj[neg_node(na)].append(nb)
        radj[nb].append(neg_node(na))
        adj[neg_node(nb)].append(na)
        radj[na].append(neg_node(nb))

    # Pass 1: DFS on original graph, record finish order
    visited = [False] * N
    finish_order = []
    for start in range(N):
        if visited[start]:
            continue
        stack = [(start, False)]
        while stack:
            node, processed = stack.pop()
            if processed:
                finish_order.append(node)
                continue
            if visited[node]:
                continue
            visited[node] = True
            stack.append((node, True))
            for nb in adj[node]:
                if not visited[nb]:
                    stack.append((nb, False))

    # Pass 2: DFS on reversed graph in reverse finish order
    comp = [-1] * N
    comp_id = 0
    for start in reversed(finish_order):
        if comp[start] != -1:
            continue
        comp_id += 1
        stack = [start]
        comp[start] = comp_id
        while stack:
            node = stack.pop()
            for nb in radj[node]:
                if comp[nb] == -1:
                    comp[nb] = comp_id
                    stack.append(nb)

    # Check satisfiability: x and ¬x in same SCC => unsatisfiable
    for i in range(n):
        if comp[2 * i] == comp[2 * i + 1]:
            return None

    # Construct assignment: x_i = TRUE if comp[x_i] > comp[¬x_i]
    # (process SCCs in reverse topological order; sink SCCs get TRUE)
    assignment = []
    for i in range(n):
        if comp[2 * i] > comp[2 * i + 1]:
            assignment.append(i + 1)
        else:
            assignment.append(-(i + 1))

    return assignment


def main():
    input_file = '2sat_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            content = f.read()
    else:
        content = sys.stdin.read()

    lines = content.strip().splitlines()
    k = int(lines[0].strip())

    results = []
    line_idx = 1
    for _ in range(k):
        # Skip blank lines
        while line_idx < len(lines) and lines[line_idx].strip() == '':
            line_idx += 1
        n, m = map(int, lines[line_idx].split())
        line_idx += 1
        clauses = []
        for _ in range(m):
            a, b = map(int, lines[line_idx].split())
            clauses.append((a, b))
            line_idx += 1

        assignment = solve_2sat(n, clauses)
        if assignment is None:
            results.append('0')
        else:
            results.append('1 ' + ' '.join(str(x) for x in assignment))

    print(' '.join(results))


if __name__ == '__main__':
    main()
