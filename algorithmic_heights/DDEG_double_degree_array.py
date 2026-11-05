# DDEG (Double-Degree Array)
# https://rosalind.info/problems/ddeg/
#
# What the algorithm is:
# The double-degree array D is defined such that D[i] is the sum of the
# degrees of all neighbors of vertex i. In other words, for each vertex,
# we look at its neighbors and add up how many neighbors each of those
# neighbors has.
#
# How it works conceptually:
# First, build an adjacency list representation of the graph and compute
# the degree of each vertex. Then, for each vertex i, iterate through its
# neighbors and sum their degrees.
#
# Time and space complexity:
# Time: O(n + m) — build adjacency list and degrees in O(n + m), then
# compute double degrees by iterating through each edge once from each
# endpoint's perspective: O(sum of degrees) = O(2m) = O(m).
# Space: O(n + m) — store adjacency list and degree array.
#
# Graph representation:
# Adjacency list using a dictionary mapping each vertex to a list of its
# neighbors.
#
# Where this appears in bioinformatics pipelines:
# Double-degree-like measures appear in network analysis of biological
# networks. For example, in protein-protein interaction networks, the
# "neighborhood connectivity" of a protein (average degree of its
# interaction partners) can indicate whether it interacts with highly
# connected hub proteins, which may be biologically significant.


def double_degree_array(n, edges):
    """Compute the double-degree array for an undirected graph.

    D[i] is the sum of the degrees of all neighbors of vertex i.

    Args:
        n: Number of vertices (numbered 1 to n).
        edges: A list of (u, v) tuples representing undirected edges.

    Returns:
        A list D where D[i] is the double degree of vertex i (D[0] unused).
    """
    # Build adjacency list and compute degrees
    adj = {i: [] for i in range(1, n + 1)}
    degree = [0] * (n + 1)

    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
        degree[u] += 1
        degree[v] += 1

    # Compute double degree for each vertex
    double_degree = [0] * (n + 1)
    for i in range(1, n + 1):
        for neighbor in adj[i]:
            double_degree[i] += degree[neighbor]

    return double_degree


def main():
    import sys
    import os

    input_file = 'ddeg_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            content = f.read()
    else:
        content = sys.stdin.read()

    lines = content.strip().splitlines()
    n, m = map(int, lines[0].split())
    edges = []
    for i in range(1, m + 1):
        u, v = map(int, lines[i].split())
        edges.append((u, v))

    dd = double_degree_array(n, edges)

    print(' '.join(str(dd[i]) for i in range(1, n + 1)))


if __name__ == '__main__':
    main()
