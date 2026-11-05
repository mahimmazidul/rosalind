# DAG (Testing Acyclicity)
# https://rosalind.info/problems/dag/
#
# What the algorithm is:
# A directed graph is acyclic (a DAG) if it contains no directed cycles.
# We can detect cycles using DFS: color vertices white (unvisited), gray
# (currently being explored), or black (fully explored). If we encounter
# a gray vertex during DFS, we've found a back edge and thus a cycle.
#
# How it works conceptually:
# Start DFS from each unvisited vertex. Mark it gray when we start
# exploring it, and black when we're done. If we find an edge to a gray
# vertex, that's a back edge → cycle exists → not a DAG.
#
# Time and space complexity:
# Time: O(k * (n + m)) — DFS visits each vertex and edge once per graph.
# Space: O(n + m) — adjacency list plus color array and recursion stack.
#
# Graph representation:
# Adjacency list using a dictionary mapping each vertex to its outgoing neighbors.
#
# Where this appears in bioinformatics:
# DAGs model many biological processes — gene regulatory networks (no
# feedback loops), metabolic pathways (reactions flow in one direction),
# phylogenetic trees (evolution flows from ancestors to descendants).
# Testing acyclicity validates these models.


import sys
import os

sys.setrecursionlimit(10000)


def is_acyclic(n, edges):
    """Test if a directed graph is acyclic using DFS coloring.

    Args:
        n: Number of vertices (numbered 1 to n).
        edges: A list of (u, v) tuples representing directed edges u -> v.

    Returns:
        True if the graph is acyclic, False otherwise.
    """
    # Build adjacency list
    adj = {i: [] for i in range(1, n + 1)}
    for u, v in edges:
        adj[u].append(v)

    # color: 0 = white (unvisited), 1 = gray (in progress), 2 = black (done)
    color = {i: 0 for i in range(1, n + 1)}

    def dfs(u):
        color[u] = 1  # gray
        for v in adj[u]:
            if color[v] == 1:
                return False  # back edge → cycle
            if color[v] == 0:
                if not dfs(v):
                    return False
        color[u] = 2  # black
        return True

    for v in range(1, n + 1):
        if color[v] == 0:
            if not dfs(v):
                return False

    return True


def main():
    input_file = 'dag_input.txt'
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

        if is_acyclic(n, edges):
            results.append('1')
        else:
            results.append('-1')

    print(' '.join(results))


if __name__ == '__main__':
    main()
