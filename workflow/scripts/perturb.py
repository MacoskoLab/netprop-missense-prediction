import pandas as pd
import numpy as np
import networkx as nx
from snakemake.script import snakemake

gene = snakemake.params["gene"]
k = snakemake.params["k"]
t = snakemake.params["t"]

network_df = pd.read_csv(snakemake.input["genie3_links"], sep="\t")

# Read AM score from file
with open(snakemake.input["am_score"]) as f:
    line = f.readline().strip().split("\t")
    am_gene, am_score = line[0], float(line[1])

assert (
    gene == am_gene
), f"Gene in config ({gene}) does not match gene in AM score file ({am_gene})"

am_scores = {gene: am_score}


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


# Workflow
G = build_graph(network_df)
Gp = apply_perturbation_graph(G, gene, am_scores, mode="knockout", k=k, t=t)
impacts = propagate_impacts_graph(Gp, {gene: 1.0}, steps=2)

# Output
out_df = pd.DataFrame(list(impacts.items()), columns=["gene", "impact"])
out_df.to_csv(snakemake.output["scores"], index=False)
