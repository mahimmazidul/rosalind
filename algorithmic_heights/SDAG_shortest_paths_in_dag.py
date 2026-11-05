# SDAG (Shortest Paths in DAG)
# https://rosalind.info/problems/sdag/
#
# What the algorithm is:
# Single-source shortest paths in a directed acyclic graph from vertex 1.
# Since a DAG has no cycles, we can topologically sort the vertices and
# relax edges in that order — each vertex is finalized before its
# outgoing edges are processed, guaranteeing correctness in O(n + m).
#
# How it works conceptually:
# 1. Topologically sort the DAG (Kahn's algorithm).
# 2. Initialize D[1] = 0, all others = infinity.
# 3. Process vertices in topological order; for each edge (u, v, w),
#    relax: D[v] = min(D[v], D[u] + w).
# 4. Unreachable vertices (D = infinity) are output as "x".
#
# Edge weights may be negative — no problem since there are no cycles.
#
# Time and space complexity:
# Time: O(n + m) — topological sort plus one pass over all edges.
# Space: O(n + m) — adjacency list plus distance array.
#
# Graph representation:
# Adjacency list of (neighbor, weight) tuples.
#
# Where this appears in bioinformatics:
# Shortest paths in DAGs model optimal sequence alignment in directed
# acyclic graphs (e.g., pangenome graphs, variant graphs), where the
# best path through a sequence graph is found in linear time.

from collections import deque
import sys
import os


def topological_sort(n, adj):
    """Topological sort of a DAG using Kahn's algorithm."""
    in_degree = [0] * (n + 1)
    for u in range(1, n + 1):
        for v, _ in adj[u]:
            in_degree[v] += 1

    queue = deque()
    for i in range(1, n + 1):
        if in_degree[i] == 0:
            queue.append(i)

    order = []
    while queue:
        u = queue.popleft()
        order.append(u)
        for v, _ in adj[u]:
            in_degree[v] -= 1
            if in_degree[v] == 0:
                queue.append(v)

    return order


def shortest_paths_dag(n, edges):
    """Shortest paths from vertex 1 in a DAG.

    Args:
        n: Number of vertices (numbered 1 to n).
        edges: A list of (u, v, w) tuples.

    Returns:
        A list D where D[i] is the shortest distance from vertex 1 to i,
        or None (rendered as 'x') if unreachable.
    """
    adj = {i: [] for i in range(1, n + 1)}
    for u, v, w in edges:
        adj[u].append((v, w))

    order = topological_sort(n, adj)

    INF = float('inf')
    dist = [INF] * (n + 1)
    dist[1] = 0

    for u in order:
        if dist[u] == INF:
            continue
        for v, w in adj[u]:
            if dist[u] + w < dist[v]:
                dist[v] = dist[u] + w

    return dist


def main():
    input_file = 'sdag_input.txt'
    if os.path.exists(input_file):
        with open(input_file, 'r') as f:
            content = f.read()
    else:
        content = sys.stdin.read()

    lines = content.strip().splitlines()
    n, m = map(int, lines[0].split())
    edges = []
    for i in range(1, m + 1):
        u, v, w = map(int, lines[i].split())
        edges.append((u, v, w))

    dist = shortest_paths_dag(n, edges)

    output = []
    for i in range(1, n + 1):
        if dist[i] == float('inf'):
            output.append('x')
        else:
            output.append(str(int(dist[i])))

    print(' '.join(output))


if __name__ == '__main__':
    main()
