# HDAG (Hamiltonian Path in DAG)
# https://rosalind.info/problems/hdag/
#
# What the algorithm is:
# A Hamiltonian path visits every vertex exactly once. In general graphs,
# finding one is NP-hard. But in a DAG, it's easy: a Hamiltonian path
# exists if and only if the topological ordering has a direct edge
# between every consecutive pair of vertices.
#
# How it works conceptually:
# Compute a topological sort of the DAG. Then check if every consecutive
# pair (v[i], v[i+1]) in the ordering has a direct edge. If all pairs
# are connected, the topological order IS a Hamiltonian path. If any
# pair is missing an edge, no Hamiltonian path exists.
#
# Why this works: In a DAG, any Hamiltonian path must respect the
# topological order (all edges go forward). So the Hamiltonian path,
# if it exists, must be exactly the topological order with all
# consecutive pairs connected by edges.
#
# Time and space complexity:
# Time: O(n + m) — topological sort plus checking n-1 consecutive pairs.
# Space: O(n + m) — adjacency list.
#
# Graph representation:
# Adjacency list plus adjacency set for O(1) edge lookups.
#
# Where this appears in bioinformatics:
# Hamiltonian paths in DAGs appear in genome assembly — finding a path
# through all contigs in an overlap graph, where each contig is visited
# exactly once.


from collections import deque
import sys
import os


def topological_sort(n, edges):
    """Return a topological ordering of a DAG using Kahn's algorithm."""
    adj = {i: [] for i in range(1, n + 1)}
    in_degree = [0] * (n + 1)
    for u, v in edges:
        adj[u].append(v)
        in_degree[v] += 1

    queue = deque()
    for i in range(1, n + 1):
        if in_degree[i] == 0:
            queue.append(i)

    result = []
    while queue:
        u = queue.popleft()
        result.append(u)
        for v in adj[u]:
            in_degree[v] -= 1
            if in_degree[v] == 0:
                queue.append(v)

    return result


def has_hamiltonian_path(n, edges):
    """Check if a DAG has a Hamiltonian path and return it if so.

    Args:
        n: Number of vertices (numbered 1 to n).
        edges: A list of (u, v) tuples representing directed edges.

    Returns:
        A tuple (found, path) where found is a boolean and path is a
        list of vertices forming the Hamiltonian path (or empty list).
    """
    # Build adjacency set for O(1) edge lookups
    adj_set = {i: set() for i in range(1, n + 1)}
    for u, v in edges:
        adj_set[u].add(v)

    order = topological_sort(n, edges)

    # Check if every consecutive pair has a direct edge
    for i in range(len(order) - 1):
        if order[i + 1] not in adj_set[order[i]]:
            return False, []

    return True, order


def main():
    input_file = 'hdag_input.txt'
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

        found, path = has_hamiltonian_path(n, edges)
        if found:
            results.append('1 ' + ' '.join(str(x) for x in path))
        else:
            results.append('-1')

    print(' '.join(results))


if __name__ == '__main__':
    main()
