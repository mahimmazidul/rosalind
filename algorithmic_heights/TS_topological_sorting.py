# TS (Topological Sorting)
# https://rosalind.info/problems/ts/
#
# What the algorithm is:
# Topological sorting orders the vertices of a directed acyclic graph
# (DAG) so that for every directed edge u → v, u comes before v in the
# ordering. It's only possible for DAGs.
#
# How it works conceptually (Kahn's algorithm):
# Repeatedly find a vertex with no incoming edges (in-degree 0), output
# it, and remove it from the graph (decrement in-degrees of its
# neighbors). Continue until all vertices are output.
#
# Think of it as course scheduling: a course can only be taken after
# all its prerequisites are completed. Output courses in an order where
# prerequisites always come first.
#
# Time and space complexity:
# Time: O(n + m) — each vertex processed once, each edge examined once.
# Space: O(n + m) — adjacency list plus in-degree array.
#
# Graph representation:
# Adjacency list using a dictionary mapping each vertex to its outgoing
# neighbors.
#
# Where this appears in bioinformatics:
# Topological sort orders steps in biological pipelines — for example,
# ordering steps in a metabolic pathway (reaction A must happen before
# reaction B), or ordering tasks in a bioinformatics workflow (quality
# control before alignment before variant calling).


from collections import deque
import sys
import os


def topological_sort(n, edges):
    """Return a topological ordering of a DAG using Kahn's algorithm.

    Args:
        n: Number of vertices (numbered 1 to n).
        edges: A list of (u, v) tuples representing directed edges
               u -> v.

    Returns:
        A list of vertices in topological order.
    """
    # Build adjacency list and in-degree count
    adj = {i: [] for i in range(1, n + 1)}
    in_degree = [0] * (n + 1)

    for u, v in edges:
        adj[u].append(v)
        in_degree[v] += 1

    # Initialize queue with all vertices of in-degree 0
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


def main():
    input_file = 'ts_input.txt'
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

    result = topological_sort(n, edges)
    print(' '.join(str(x) for x in result))


if __name__ == '__main__':
    main()
