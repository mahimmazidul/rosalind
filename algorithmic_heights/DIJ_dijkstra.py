# DIJ (Dijkstra's Algorithm)
# https://rosalind.info/problems/dij/
#
# What the algorithm is:
# Dijkstra's algorithm computes the shortest path from a source vertex to
# all other vertices in a weighted graph with non-negative edge weights.
# It greedily selects the closest unvisited vertex and relaxes its edges.
#
# How it works conceptually:
# Maintain a priority queue of (distance, vertex) pairs. Start with the
# source at distance 0. Repeatedly extract the closest vertex, and for
# each neighbor, check if going through the current vertex gives a shorter
# path. If so, update the distance and add to the queue.
#
# Time and space complexity:
# Time: O((n + m) log n) — each edge may cause a heap push, each vertex
# extracted once.
# Space: O(n + m) — adjacency list plus distance array and heap.
#
# Graph representation:
# Adjacency list using a dictionary mapping each vertex to a list of
# (neighbor, weight) tuples.
#
# Where this appears in bioinformatics:
# Dijkstra's finds shortest paths in weighted biological networks — for
# example, minimum-cost paths in sequence alignment graphs, or shortest
# paths in phylogenetic trees with branch lengths.


import heapq
import sys
import os


def dijkstra(n, edges):
    """Compute shortest distances from vertex 1 in a directed weighted graph.

    Args:
        n: Number of vertices (numbered 1 to n).
        edges: A list of (u, v, w) tuples representing directed edges
               u -> v with weight w.

    Returns:
        A list D where D[i] is the shortest distance from vertex 1 to
        vertex i, or -1 if unreachable. D[0] is unused.
    """
    # Build adjacency list
    adj = {i: [] for i in range(1, n + 1)}
    for u, v, w in edges:
        adj[u].append((v, w))

    # Dijkstra
    INF = float('inf')
    distance = [INF] * (n + 1)
    distance[1] = 0
    heap = [(0, 1)]

    while heap:
        d, u = heapq.heappop(heap)
        if d > distance[u]:
            continue  # stale entry
        for v, w in adj[u]:
            new_dist = d + w
            if new_dist < distance[v]:
                distance[v] = new_dist
                heapq.heappush(heap, (new_dist, v))

    # Convert INF to -1
    for i in range(1, n + 1):
        if distance[i] == INF:
            distance[i] = -1

    return distance


def main():
    input_file = 'dij_input.txt'
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

    distances = dijkstra(n, edges)

    print(' '.join(str(distances[i]) for i in range(1, n + 1)))


if __name__ == '__main__':
    main()
