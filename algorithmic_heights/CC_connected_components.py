# CC (Connected Components)
# https://rosalind.info/problems/cc/
#
# What the algorithm is:
# Connected components are maximal subgraphs where every pair of vertices
# is connected by a path. To count them, use DFS: start from an unvisited
# vertex, explore all reachable vertices (marking them visited), that's one
# component. Repeat until all vertices are visited.
#
# How it works conceptually:
# Think of the graph as islands. Each island is a connected component.
# Pick an unvisited vertex, walk along edges to visit everything reachable
# from it (that's one island), mark them all visited, then pick another
# unvisited vertex and repeat. The number of times you start a new walk
# equals the number of connected components.
#
# Time and space complexity:
# Time: O(n + m) — each vertex visited once, each edge examined once.
# Space: O(n + m) — adjacency list plus visited set and recursion stack.
#
# Graph representation:
# Adjacency list using a dictionary.
#
# Where this appears in bioinformatics pipelines:
# Connected components appear in biological network analysis — for
# example, finding clusters of interacting proteins, or grouping genes
# that are co-expressed across multiple conditions.


import sys
import os

sys.setrecursionlimit(10000)


def count_connected_components(n, edges):
    """Count connected components in an undirected graph using DFS.

    Args:
        n: Number of vertices (numbered 1 to n).
        edges: A list of (u, v) tuples representing undirected edges.

    Returns:
        The number of connected components.
    """
    # Build adjacency list
    adj = {i: [] for i in range(1, n + 1)}
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)

    visited = set()
    components = 0

    def dfs(v):
        visited.add(v)
        for neighbor in adj[v]:
            if neighbor not in visited:
                dfs(neighbor)

    for v in range(1, n + 1):
        if v not in visited:
            components += 1
            dfs(v)

    return components


def main():
    input_file = 'cc_input.txt'
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

    result = count_connected_components(n, edges)
    print(result)


if __name__ == '__main__':
    main()
