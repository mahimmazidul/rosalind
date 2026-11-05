# BIP (Testing Bipartiteness)
# https://rosalind.info/problems/bip/
#
# What the algorithm is:
# A graph is bipartite if its vertices can be divided into two sets such
# that no edge connects vertices within the same set. Equivalently, the
# graph can be 2-colored: assign each vertex one of two colors so that
# adjacent vertices always have different colors.
#
# How it works conceptually:
# Use BFS or DFS to 2-color the graph. Start from an uncolored vertex,
# assign it color 0, then assign all its neighbors color 1, then their
# neighbors color 0, and so on. If we ever find an edge connecting two
# vertices of the same color, the graph is not bipartite.
#
# Time and space complexity:
# Time: O(k * (n + m)) — for each graph, BFS visits each vertex and edge once.
# Space: O(n + m) — adjacency list plus color array.
#
# Graph representation:
# Adjacency list using a dictionary.
#
# Where this appears in bioinformatics pipelines:
# Bipartite graphs model relationships between two different types of
# entities — for example, genes and diseases (a gene is associated with
# a disease, but genes aren't associated with other genes, and diseases
# aren't associated with other diseases). Testing bipartiteness validates
# such data models.


import sys
import os

sys.setrecursionlimit(10000)


def is_bipartite(n, edges):
    """Test if an undirected graph is bipartite using BFS 2-coloring.

    Args:
        n: Number of vertices (numbered 1 to n).
        edges: A list of (u, v) tuples representing undirected edges.

    Returns:
        True if the graph is bipartite, False otherwise.
    """
    # Build adjacency list
    adj = {i: [] for i in range(1, n + 1)}
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)

    # color[v] = 0 or 1, or None if uncolored
    color = {}

    from collections import deque

    for start in range(1, n + 1):
        if start in color:
            continue
        # Start BFS from this uncolored vertex
        color[start] = 0
        queue = deque([start])

        while queue:
            u = queue.popleft()
            for v in adj[u]:
                if v not in color:
                    color[v] = 1 - color[u]
                    queue.append(v)
                elif color[v] == color[u]:
                    return False

    return True


def main():
    input_file = 'bip_input.txt'
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

        if is_bipartite(n, edges):
            results.append('1')
        else:
            results.append('-1')

    print(' '.join(results))


if __name__ == '__main__':
    main()
