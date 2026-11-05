# SC (Semi-Connected Graph)
# https://rosalind.info/problems/sc/
#
# What the algorithm is:
# Check if a directed graph is semi-connected: for every pair of vertices
# (u, v), either u reaches v or v reaches u.
#
# How it works conceptually:
# 1. Find all strongly connected components (SCCs) using Kosaraju's algo.
#    Within an SCC, all vertices are mutually reachable.
# 2. Condense the graph into a DAG where each node is an SCC.
# 3. The graph is semi-connected iff the condensed DAG is "chain-like":
#    there exists a topological ordering where every consecutive pair of
#    SCCs has a direct edge (i.e., the condensed DAG has a Hamiltonian
#    path). This guarantees every pair of SCCs is comparable by reachability.
# 4. Topological sort the condensed DAG, then verify every consecutive
#    pair has an edge.
#
# Time and space complexity:
# Time: O(n + m) — SCC finding plus condensation plus topological check.
# Space: O(n + m) — adjacency lists.
#
# Graph representation:
# Adjacency lists for original graph, reverse graph, and condensed DAG.
#
# Where this appears in bioinformatics:
# Semi-connectivity checks appear in hierarchical pathway analysis —
# verifying that a biological network has a total order of influence
# (every module can affect or be affected by every other).

import sys
import os


def kosaraju_scc(n, edges):
    """Find SCCs using Kosaraju's algorithm. Returns comp array and number of SCCs."""
    adj = {i: [] for i in range(1, n + 1)}
    radj = {i: [] for i in range(1, n + 1)}
    for u, v in edges:
        adj[u].append(v)
        radj[v].append(u)

    # Pass 1: finish order on original graph
    visited = [False] * (n + 1)
    finish_order = []
    for start in range(1, n + 1):
        if visited[start]:
            continue
        stack = [(start, False)]
        while stack:
            node, processed = stack.pop()
            if processed:
                finish_order.append(node)
                continue
            if visited[node]:
                continue
            visited[node] = True
            stack.append((node, True))
            for nb in adj[node]:
                if not visited[nb]:
                    stack.append((nb, False))

    # Pass 2: assign components on reversed graph
    comp = [-1] * (n + 1)
    comp_id = 0
    for start in reversed(finish_order):
        if comp[start] != -1:
            continue
        comp_id += 1
        stack = [start]
        comp[start] = comp_id
        while stack:
            node = stack.pop()
            for nb in radj[node]:
                if comp[nb] == -1:
                    comp[nb] = comp_id
                    stack.append(nb)

    return comp, comp_id


def is_semi_connected(n, edges):
    """Check if a directed graph is semi-connected.

    Args:
        n: Number of vertices (numbered 1 to n).
        edges: A list of (u, v) tuples representing directed edges.

    Returns:
        True if semi-connected, False otherwise.
    """
    comp, num_scc = kosaraju_scc(n, edges)

    # Build condensed DAG
    cadj = {i: set() for i in range(1, num_scc + 1)}
    for u, v in edges:
        cu, cv = comp[u], comp[v]
        if cu != cv:
            cadj[cu].add(cv)

    # Topological sort of condensed DAG (Kahn's algorithm)
    indeg = {i: 0 for i in range(1, num_scc + 1)}
    for u in range(1, num_scc + 1):
        for v in cadj[u]:
            indeg[v] += 1

    from collections import deque
    queue = deque()
    for i in range(1, num_scc + 1):
        if indeg[i] == 0:
            queue.append(i)

    topo = []
    while queue:
        u = queue.popleft()
        topo.append(u)
        for v in cadj[u]:
            indeg[v] -= 1
            if indeg[v] == 0:
                queue.append(v)

    # Check every consecutive pair has a direct edge
    for i in range(len(topo) - 1):
        if topo[i + 1] not in cadj[topo[i]]:
            return False

    return True


def main():
    input_file = 'sc_input.txt'
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

        if is_semi_connected(n, edges):
            results.append('1')
        else:
            results.append('-1')

    print(' '.join(results))


if __name__ == '__main__':
    main()
