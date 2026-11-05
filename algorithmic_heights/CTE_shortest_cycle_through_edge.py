# CTE (Shortest Cycle Through Edge)
# https://rosalind.info/problems/cte/
#
# What the algorithm is:
# A cycle through edge (u, v) consists of the edge itself plus a path
# from v back to u. The shortest cycle through (u, v) is therefore
# weight(u,v) + shortest_path(v, u). If no path exists from v to u,
# there's no cycle through this edge.
#
# How it works conceptually:
# For each graph, take the first edge (u, v, w). Temporarily remove it
# from the graph. Run Dijkstra's algorithm from v to find the shortest
# path to u. If a path exists, the answer is w + dist(v, u). Otherwise,
# output -1.
#
# Time and space complexity:
# Time: O(k * (m log n)) — Dijkstra for each graph.
# Space: O(n + m) — adjacency list per graph.
#
# Graph representation:
# Adjacency list using a dictionary mapping each vertex to a list of
# (neighbor, weight) tuples.
#
# Where this appears in bioinformatics:
# Finding shortest cycles through specific edges appears in network
# reliability analysis — for example, in gene regulatory networks,
# finding alternative feedback paths when a regulatory edge is disrupted.


import heapq
import sys
import os


def shortest_path_from_to(n, adj, source, target):
    """Dijkstra: shortest distance from source to target.

    Args:
        n: Number of vertices.
        adj: Adjacency list {vertex: [(neighbor, weight), ...]}.
        source: Starting vertex.
        target: Target vertex.

    Returns:
        Shortest distance, or None if unreachable.
    """
    INF = float('inf')
    distance = [INF] * (n + 1)
    distance[source] = 0
    heap = [(0, source)]

    while heap:
        d, u = heapq.heappop(heap)
        if d > distance[u]:
            continue
        if u == target:
            return d
        for v, w in adj[u]:
            new_dist = d + w
            if new_dist < distance[v]:
                distance[v] = new_dist
                heapq.heappush(heap, (new_dist, v))

    return None if distance[target] == INF else distance[target]


def shortest_cycle_through_first_edge(n, edges):
    """Find the shortest cycle through the first edge.

    Args:
        n: Number of vertices.
        edges: A list of (u, v, w) tuples representing directed edges.

    Returns:
        The length of the shortest cycle through the first edge,
        or -1 if no such cycle exists.
    """
    if not edges:
        return -1

    u, v, w = edges[0]

    # Build adjacency list without the first edge
    adj = {i: [] for i in range(1, n + 1)}
    for eu, ev, ew in edges[1:]:
        adj[eu].append((ev, ew))

    # Find shortest path from v back to u
    dist = shortest_path_from_to(n, adj, v, u)

    if dist is None:
        return -1
    return w + dist


def main():
    input_file = 'cte_input.txt'
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
            u, v, w = map(int, lines[line_idx].split())
            edges.append((u, v, w))
            line_idx += 1

        result = shortest_cycle_through_first_edge(n, edges)
        results.append(str(result))

    print(' '.join(results))


if __name__ == '__main__':
    main()
