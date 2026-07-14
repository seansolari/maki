
from cdbg import ColoredDeBruijnGraph


def plot_cdbg(cdbg, outfile, title="cDBG", figsize=(8, 6)):
    import matplotlib.pyplot as plt
    import matplotlib.lines as pll
    import networkx as nx
    from matplotlib.colors import to_hex

    G = cdbg.graph
    pos = nx.spring_layout(G, seed=41)

    # Collect all colours
    all_colors = set()
    for _, _, data in G.edges(data=True):
        all_colors.add(data["color"])

    cmap = plt.get_cmap("tab10")
    color_list = list(all_colors)
    color_map = {c: to_hex(cmap(i % 10)) for i, c in enumerate(color_list)}

    plt.figure(figsize=figsize)

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, node_size=400, node_color="#eeeeee", edgecolors="black")
    nx.draw_networkx_labels(G, pos, font_size=8, font_family="monospace")

    # Draw edges per colour
    for u, v, key, data in G.edges(keys=True, data=True):
        color = data["color"]
        cov = cdbg.get_edge_coverage(u, v, color)

        nx.draw_networkx_edges(
            G,
            pos,
            edgelist=[(u, v)],
            connectionstyle="arc3,rad=0.2",  # separates parallel edges
            edge_color=color_map[color],
            width=1 + cov * 0.5,
            alpha=0.9,
            arrows=True
        )

    # Legend
    handles = [
        pll.Line2D([0], [0], color=color_map[c], lw=3, label=c)
        for c in color_list
    ]
    plt.legend(handles=handles)

    plt.title(title)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(outfile)
    

if __name__ == "__main__":
    cdbg = ColoredDeBruijnGraph(k=4)

    genome1 = "ATGCGATGAC"
    genome2 = "ATGCGTTGAC"

    cdbg.add_sequence(genome1, "g1")
    cdbg.add_sequence(genome2, "g2")

    plot_cdbg(cdbg, "test.pdf", title="Compacted cDBG")
    