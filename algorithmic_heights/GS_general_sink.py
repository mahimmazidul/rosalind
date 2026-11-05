# GS (General Sink / Mother Vertex)
# https://rosalind.info/problems/gs/
#
# What the algorithm is:
# Find a "mother vertex" in a directed graph — a vertex from which all
# other vertices are reachable via directed paths.
#
# How it works conceptually:
# 1. Run DFS on the graph, tracking the last vertex to finish.
# 2. If a mother vertex exists, the last finished vertex must be one
#    (intuition: DFS starting from a mother vertex reaches everything,
#    so it finishes last; if we started elsewhere and still reached the
#    mother vertex last, it's still a valid candidate).
# 3. Verify by running DFS/BFS from the candidate and checking that all
#    vertices are visited. If not, no mother vertex exists.
#
# Time and space complexity:
# Time: O(n + m) — two DFS traversals.
# Space: O(n + m) — adjacency list plus visited array.
#
# Graph representation:
# Adjacency list.
#
# Where this appears in bioinformatics:
# Mother vertex / reachability analysis appears in signaling pathway
# analysis — finding a master regulator (e.g. a transcription factor)
# from which all genes in a regulatory network are reachable.

import sys
import os


def dfs_reachable(start, adj, n):
    """Return the set of vertices reachable from start via directed edges."""
    visited = [False] * (n + 1)
    stack = [start]
    visited[start] = True
    while stack:
        node = stack.pop()
        for nb in adj[node]:
            if not visited[nb]:
                visited[nb] = True
                stack.append(nb)
    return visited


def find_mother_vertex(n, edges):
    """Find a vertex from which all others are reachable, or -1 if none.

    Args:
        n: Number of vertices (numbered 1 to n).
        edges: A list of (u, v) tuples representing directed edges.

    Returns:
        A mother vertex (int) or -1 if none exists.
    """
    adj = {i: [] for i in range(1, n + 1)}
    for u, v in edges:
        adj[u].append(v)

    # Step 1: DFS to find the last finished vertex (candidate mother vertex)
    visited = [False] * (n + 1)
    last_finished = -1

    for start in range(1, n + 1):
        if visited[start]:
            continue
        stack = [(start, False)]
        while stack:
            node, processed = stack.pop()
            if processed:
                last_finished = node
                continue
            if visited[node]:
                continue
            visited[node] = True
            stack.append((node, True))
            for nb in adj[node]:
                if not visited[nb]:
                    stack.append((nb, False))

    # Step 2: Verify the candidate reaches all vertices
    if last_finished == -1:
        return -1

    reachable = dfs_reachable(last_finished, adj, n)
    if all(reachable[1:]):
        return last_finished
    return -1


def main():
    input_file = 'gs_input.txt'
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

        mv = find_mother_vertex(n, edges)
        results.append(str(mv))

    print(' '.join(results))


if __name__ == '__main__':
    main()
