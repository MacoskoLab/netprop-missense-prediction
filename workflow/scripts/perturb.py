import pandas as pd
import numpy as np
import networkx as nx

# Sample input: GENIE3 output table
data = {
    "regulatoryGene": ["Gene2", "Gene2", "Gene2", "Gene7", "Gene7"],
    "targetGene": ["Gene14", "Gene7", "Gene12", "Gene18", "Gene17"],
    "weight": [0.0, 0.0, 0.0, 0.7820, 0.6928],
}
network_df = pd.DataFrame(data)

# Sample AlphaMissense pathogenicity scores
am_scores = {
    "Gene2": 0.9,
    "Gene7": 0.4,
    "Gene14": 0.1,
    "Gene12": 0.2,
    "Gene18": 0.3,
    "Gene17": 0.5,
}


def logistic_transform(score, k=10, t=0.5):
    return 1 / (1 + np.exp(-k * (score - t)))


def build_graph(df):
    G = nx.DiGraph()
    for _, row in df.iterrows():
        G.add_edge(row["regulatoryGene"], row["targetGene"], weight=row["weight"])
    return G


def apply_perturbation_graph(G, gene, am_scores, mode="knockout", k=10, t=0.5):
    H = G.copy()
    severity = logistic_transform(am_scores.get(gene, 0), k, t)
    if mode == "knockout":
        H.remove_edges_from(list(H.out_edges(gene)))
    elif mode == "knockdown":
        scale = 1 - severity
        for u, v, data in H.out_edges(gene, data=True):
            data["weight"] *= scale
    else:
        raise ValueError("mode must be 'knockout' or 'knockdown'")
    return H


def propagate_impacts_graph(G, initial_impacts, steps=2):
    nodes = list(G.nodes)
    idx = {g: i for i, g in enumerate(nodes)}
    adj = nx.to_numpy_array(G, nodelist=nodes, weight="weight")
    vec = np.zeros(len(nodes))
    for g, val in initial_impacts.items():
        vec[idx[g]] = val
    for _ in range(steps):
        vec = vec @ adj
    return {nodes[i]: vec[i] for i in range(len(nodes))}


# Example workflow
G = build_graph(network_df)
Gp = apply_perturbation_graph(G, "Gene2", am_scores, mode="knockout")
impacts = propagate_impacts_graph(Gp, {"Gene2": 1.0}, steps=2)

# Output
print("Nodes:", Gp.nodes)
print("Edges with weights:", list(Gp.edges(data=True)))
print("Propagated impacts:", impacts)
