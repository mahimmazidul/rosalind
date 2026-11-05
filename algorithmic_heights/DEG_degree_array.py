# DEG (Degree Array)
# https://rosalind.info/problems/deg/
#
# What the algorithm is:
# In an undirected graph, the degree of a vertex is the number of edges
# incident to it — equivalently, the number of neighbors it has. The degree
# array lists the degree of every vertex in the graph.
#
# How it works conceptually:
# Read the graph as an edge list. For each edge (u, v), increment the degree
# counter of both u and v (since the graph is undirected, the edge contributes
# to the degree of both endpoints). After processing all edges, the degree
# array contains the answer.
#
# Time and space complexity:
# Time: O(n + m) — initialize n counters to zero, then process m edges.
# Space: O(n) — store one integer per vertex.
#
# Graph representation:
# The graph is given as an edge list. We use an array of size n+1 (1-indexed)
# to store degrees, since vertices are numbered 1 through n.
#
# Where this appears in bioinformatics pipelines:
# Graphs model many biological networks — protein-protein interaction
# networks, gene regulatory networks, metabolic pathways. The degree of a
# node in a PPI network, for example, indicates how many other proteins it
# interacts with, which can suggest its importance (hub proteins often have
# high degree and are essential for cell survival).


def degree_array(n, edges):
    """Compute the degree of each vertex in an undirected graph.

    Args:
        n: Number of vertices (numbered 1 to n).
        edges: A list of (u, v) tuples representing undirected edges.

    Returns:
        A list D where D[i] is the degree of vertex i (D[0] is unused).
    """
    degree = [0] * (n + 1)  # 1-indexed, degree[0] unused

    for u, v in edges:
        degree[u] += 1
        degree[v] += 1

    return degree


def main():
    import sys
    import os

    input_file = 'deg_input.txt'
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

    degrees = degree_array(n, edges)

    # Output degrees of vertices 1 through n
    print(' '.join(str(degrees[i]) for i in range(1, n + 1)))


if __name__ == '__main__':
    main()
