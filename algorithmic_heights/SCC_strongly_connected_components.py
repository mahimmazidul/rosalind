# SCC (Strongly Connected Components)
# https://rosalind.info/problems/scc/
#
# What the algorithm is:
# Kosaraju's algorithm for counting strongly connected components (SCCs)
# in a directed graph. Two vertices are strongly connected if each is
# reachable from the other.
#
# How it works conceptually:
# Pass 1: DFS on the original graph, pushing vertices onto a stack in
# order of finishing time (post-order).
# Pass 2: Pop vertices from the stack and run DFS on the reversed graph
# (all edges flipped). Each DFS tree discovered in pass 2 is one SCC.
#
# Why this works: The vertex that finishes last in pass 1 belongs to a
# source SCC in the condensed graph. Processing in reverse finish order
# on the reversed graph explores exactly one SCC per DFS tree.
#
# Time and space complexity:
# Time: O(n + m) — two DFS traversals.
# Space: O(n + m) — adjacency lists for graph and reversed graph.
#
# Graph representation:
# Adjacency list for the original graph and the reversed graph.
#
# Where this appears in bioinformatics:
# SCC decomposition is used in analyzing regulatory networks — finding
# feedback loops in gene regulatory or protein interaction networks,
# where a strongly connected component represents a set of molecules
# that mutually regulate each other.

import sys
import os


def count_scc(n, edges):
    """Count the number of strongly connected components using Kosaraju's algorithm.

    Args:
        n: Number of vertices (numbered 1 to n).
        edges: A list of (u, v) tuples representing directed edges.

    Returns:
        The number of strongly connected components.
    """
    # Build adjacency lists for original and reversed graphs
    adj = {i: [] for i in range(1, n + 1)}
    radj = {i: [] for i in range(1, n + 1)}
    for u, v in edges:
        adj[u].append(v)
        radj[v].append(u)

    # Pass 1: DFS on original graph, record finish order
    visited = [False] * (n + 1)
    finish_order = []

    for start in range(1, n + 1):
        if visited[start]:
            continue
        # Iterative DFS to avoid recursion limit
        stack = [(start, False)]
        while stack:
            node, processed = stack.pop()
            if processed:
                finish_order.append(node)
                continue
            if visited[node]:
                continue
            visited[node] = True
            stack.append((node, True))  # will be added to finish_order after children
            for nb in adj[node]:
                if not visited[nb]:
                    stack.append((nb, False))

    # Pass 2: DFS on reversed graph in reverse finish order
    visited = [False] * (n + 1)
    scc_count = 0

    for start in reversed(finish_order):
        if visited[start]:
            continue
        scc_count += 1
        stack = [start]
        visited[start] = True
        while stack:
            node = stack.pop()
            for nb in radj[node]:
                if not visited[nb]:
                    visited[nb] = True
                    stack.append(nb)

    return scc_count


def main():
    input_file = 'scc_input.txt'
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

    print(count_scc(n, edges))


if __name__ == '__main__':
    main()
