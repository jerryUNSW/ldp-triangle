import argparse
import networkx as nx
import random

def generate_graph(n, m):
    if m > (n * (n - 1)) // 2:
        raise ValueError("Too many edges for the number of vertices")

    G = nx.Graph()

    # Add nodes
    G.add_nodes_from(range(n))

    # Add edges randomly
    edges = set()
    while len(edges) < m:
        u = random.randint(0, n - 1)
        v = random.randint(0, n - 1)
        if u != v and (u, v) not in edges and (v, u) not in edges:
            G.add_edge(u, v)
            edges.add((u, v))

    return G

def save_graph(G, filename):
    n = G.number_of_nodes()
    m = G.number_of_edges()
    with open(filename, 'w') as f:
        f.write(f"{n} {m}\n")
        for u, v in G.edges():
            f.write(f"{u} {v}\n")

def main():
    parser = argparse.ArgumentParser(description="Generate an undirected, unweighted graph dataset.")
    parser.add_argument('n', type=int, help="Number of vertices")
    parser.add_argument('m', type=int, help="Number of edges")
    args = parser.parse_args()

    graph = generate_graph(args.n, args.m)
    save_graph(graph, "test2.edges")

if __name__ == "__main__":
    main()
