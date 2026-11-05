# BF (Bellman-Ford Algorithm)
# https://rosalind.info/problems/bf/
#
# What the algorithm is:
# Bellman-Ford computes shortest paths from a single source in a weighted
# directed graph, even when edges have negative weights. It works by
# relaxing all edges repeatedly — after n-1 rounds, all shortest paths
# are found (assuming no negative cycles).
#
# How it works conceptually:
# Initialize all distances to infinity except the source (distance 0).
# Then repeat n-1 times: for each edge (u, v, w), if distance[u] + w <
# distance[v], update distance[v]. After n-1 rounds, no more improvements
# are possible (unless there's a negative cycle).
#
# Time and space complexity:
# Time: O(n * m) — n-1 rounds, each examining all m edges.
# Space: O(n + m) — distance array plus edge list.
#
# Where this appears in bioinformatics:
# Bellman-Ford handles graphs with negative weights, which appear in
# sequence alignment scoring (gap penalties can be negative) and in
# some phylogenetic network analyses where edge weights represent
# evolutionary distances that can be negative.


import sys
import os


def bellman_ford(n, edges):
    """Compute shortest distances from vertex 1 using Bellman-Ford.

    Args:
        n: Number of vertices (numbered 1 to n).
        edges: A list of (u, v, w) tuples representing directed edges
               u -> v with weight w.

    Returns:
        A list D where D[i] is the shortest distance from vertex 1 to
        vertex i, or None if unreachable. D[0] is unused.
    """
    INF = float('inf')
    distance = [INF] * (n + 1)
    distance[1] = 0

    # Relax all edges n-1 times
    for _ in range(n - 1):
        updated = False
        for u, v, w in edges:
            if distance[u] != INF and distance[u] + w < distance[v]:
                distance[v] = distance[u] + w
                updated = True
        if not updated:
            break  # Early exit if no changes

    return distance


def main():
    input_file = 'bf_input.txt'
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

    distances = bellman_ford(n, edges)

    # Format output: use 'x' for unreachable
    output = []
    for i in range(1, n + 1):
        if distances[i] == float('inf'):
            output.append('x')
        else:
            output.append(str(distances[i]))

    print(' '.join(output))


if __name__ == '__main__':
    main()
