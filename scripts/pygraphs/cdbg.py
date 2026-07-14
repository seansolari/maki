
from collections import defaultdict
import networkx as nx


class ColoredDeBruijnGraph:
    def __init__(self, k: int):
        self.k = k
        self.graph = nx.MultiDiGraph()

        # Node colours
        self.node_colors = defaultdict(set)

        # Edge coverage: (u, v, color) -> count
        self.edge_coverage = defaultdict(int)

    def _kmer_to_edge(self, kmer):
        return kmer[:-1], kmer[1:]

    def add_sequence(self, sequence: str, color: str):
        for i in range(len(sequence) - self.k + 1):
            self._add_kmer(sequence[i:i + self.k], color)

    def add_read(self, read: str, color: str):
        for i in range(len(read) - self.k + 1):
            self._add_kmer(read[i:i + self.k], color)

    def _add_kmer(self, kmer: str, color: str):
        u, v = self._kmer_to_edge(kmer)

        # Add nodes
        self.graph.add_node(u)
        self.graph.add_node(v)

        self.node_colors[u].add(color)
        self.node_colors[v].add(color)

        edge_key = color  # one edge per color

        if self.graph.has_edge(u, v, key=edge_key):
            self.edge_coverage[(u, v, color)] += 1
        else:
            self.graph.add_edge(u, v, key=edge_key, color=color)
            self.edge_coverage[(u, v, color)] = 1

    def get_edge_coverage(self, u, v, color):
        return self.edge_coverage.get((u, v, color), 0)

    def get_colors_between(self, u, v):
        """Return all colours with edges between u and v."""
        if self.graph.has_edge(u, v):
            return list(self.graph[u][v].keys())  # keys = colours
        return []

