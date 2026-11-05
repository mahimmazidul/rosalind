# BFS (Breadth-First Search)
# https://rosalind.info/problems/bfs/
#
# What the algorithm is:
# Breadth-First Search (BFS) explores a graph level by level, starting from
# a source vertex. It visits all vertices at distance 1, then all at
# distance 2, and so on. For unweighted graphs, BFS finds the shortest
# path from the source to every reachable vertex.
#
# How it works conceptually:
# Use a queue. Start with the source vertex at distance 0. Dequeue a
# vertex, examine all its outgoing neighbors. Any unvisited neighbor gets
# distance = current distance + 1 and is enqueued. Repeat until the queue
# is empty. Vertices never reached stay at distance -1.
#
# Time and space complexity:
# Time: O(n + m) — each vertex is enqueued once, each edge examined once.
# Space: O(n + m) — adjacency list plus distance and queue storage.
#
# Graph representation:
# Adjacency list using a dictionary mapping each vertex to a list of
# outgoing neighbors.
#
# Where this appears in bioinformatics pipelines:
# BFS finds shortest paths in unweighted biological networks — for example,
# the minimum number of interactions between two proteins in a PPI network,
# or the fewest genetic changes between two genotypes in a mutation graph.


from collections import deque
import sys
import os


def bfs_shortest_distances(n, edges):
    """Compute shortest distances from vertex 1 in a directed graph.

    Args:
        n: Number of vertices (numbered 1 to n).
        edges: A list of (u, v) tuples representing directed edges u -> v.

    Returns:
        A list D where D[i] is the shortest distance from vertex 1 to
        vertex i, or -1 if unreachable. D[0] is unused.
    """
    # Build adjacency list
    adj = {i: [] for i in range(1, n + 1)}
    for u, v in edges:
        adj[u].append(v)

    # BFS
    distance = [-1] * (n + 1)
    distance[1] = 0
    queue = deque([1])

    while queue:
        u = queue.popleft()
        for v in adj[u]:
            if distance[v] == -1:
                distance[v] = distance[u] + 1
                queue.append(v)

    return distance


def main():
    input_file = 'bfs_input.txt'
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

    distances = bfs_shortest_distances(n, edges)

    print(' '.join(str(distances[i]) for i in range(1, n + 1)))


if __name__ == '__main__':
    main()
