import sys
import os

# SQ (Square in Graph)
# https://rosalind.info/problems/sq/
#
# What the algorithm is:
# A graph contains a simple cycle of length 4 if and only if there exist
# two distinct vertices that share at least two common neighbors. If
# vertices u and v both connect to a and b, then u-a-v-b-u forms a
# 4-cycle.
#
# How it works conceptually:
# For each vertex u, look at all pairs of its neighbors (a, b). If the
# same pair (a, b) is found as neighbors of two different vertices u and
# v, then u-a-v-b-u is a 4-cycle. We use a dictionary to map each
# neighbor pair to the vertex where it was first seen.
#
# Time and space complexity:
# Time: O(sum of degree^2) — for each vertex, examine all pairs of its
# neighbors. Worst case O(n * d^2) where d is max degree.
# Space: O(n * d^2) — store neighbor pairs in a dictionary.
#
# Graph representation:
# Adjacency list using a dictionary mapping each vertex to a set of
# neighbors.
#
# Where this appears in bioinformatics:
# Finding 4-cycles (or other small cycles) in biological networks can
# reveal functional modules — for example, in protein-protein interaction
# networks, small cycles might indicate protein complexes.


def has_four_cycle(n, edges):
    """Test if an undirected graph contains a simple 4-cycle.

    Args:
        n: Number of vertices (numbered 1 to n).
        edges: A list of (u, v) tuples representing undirected edges.

    Returns:
        True if the graph contains a 4-cycle, False otherwise.
    """
    # Build adjacency list as sets
    adj = {i: set() for i in range(1, n + 1)}
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)

    # Map each neighbor pair (a, b) to the vertex where it was first seen
    pair_seen = {}

    for u in range(1, n + 1):
        neighbors = sorted(adj[u])
        # Check all pairs of neighbors of u
        for i in range(len(neighbors)):
            for j in range(i + 1, len(neighbors)):
                a = neighbors[i]
                b = neighbors[j]
                pair = (a, b)
                if pair in pair_seen:
                    # Same pair seen at a different vertex → 4-cycle
                    if pair_seen[pair] != u:
                        return True
                else:
                    pair_seen[pair] = u

    return False


def main():
    input_file = 'sq_input.txt'
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
        edges = []
        for _ in range(m):
            u, v = map(int, lines[line_idx].split())
            edges.append((u, v))
            line_idx += 1

        if has_four_cycle(n, edges):
            results.append('1')
        else:
            results.append('-1')

    print(' '.join(results))


if __name__ == '__main__':
    main()
