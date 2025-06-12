from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.script import Snakemake

    snakemake: Snakemake
    snakemake = None  # type: ignore

import pandas as pd
import numpy as np
import networkx as nx

# Snakemake parameters
steps = snakemake.params.get("steps", 2)
k = snakemake.params.get("k", 10.0)
t = snakemake.params.get("t", 0.5)

# Select pathogenicity score transform: "threshold" or "sigmoid"
transform_method = snakemake.params.get("pathogenicity_score_transform", "threshold")

# Snakemake inputs
genie3_links_file = snakemake.input["genie3_links"]
am_scores_file = snakemake.input["am_scores"]

# Snakemake outputs
scores_out_file = snakemake.output["scores"]

# Read GENIE3 links
network_df = pd.read_csv(genie3_links_file, sep="\t")

# Read AM scores: expects columns sample, gene, variant, score
am_df = pd.read_csv(am_scores_file, sep="\t", header=None)
am_df.columns = ["sample", "gene", "variant", "score"]


def threshold_grn_weights(
    W: pd.Series,
    S: pd.Series,
    threshold: float = 0.6,
) -> pd.Series:
    """Abdo 'intervention' method for GRN weights based on pathogenicity scores.
    Doing this element-wise is slow but I will come back to optimize it later.
    """
    # align indices
    W_aligned, S_aligned = W.align(S, join="inner")
    sign = np.sign(W_aligned.to_numpy())
    mag = np.abs(W_aligned.to_numpy())
    mag_new = mag.copy()
    # apply threshold per-score
    for i, score in enumerate(S_aligned.values):
        if not np.isnan(score) and score >= threshold:
            frac = (score - threshold) / (1.0 - threshold)
            frac = min(max(frac, 0.0), 1.0)
            mag_new[i] = mag[i] * (1.0 - frac)
    return pd.Series(data=sign * mag_new, index=W_aligned.index, name="W_thresholded")


def logistic_transform(score, k=10.0, t=0.5):
    """Apply logistic transformation to a score."""
    return 1 / (1 + np.exp(-k * (score - t)))


def sigmoid_grn_weights(
    W: pd.Series,
    S: pd.Series,
    k: float = k,
    t: float = t,
) -> pd.Series:
    """Sigmoid-based method for GRN weights based on pathogenicity scores.
    This method scales the weights by (1 - severity) where severity is computed
    from the logistic transformation of the score."""
    # align indices
    W_aligned, S_aligned = W.align(S, join="inner")
    # compute severity per score
    severities = S_aligned.apply(lambda score: logistic_transform(score, k, t))
    # scale weights by (1 - severity)
    W_new = W_aligned * (1.0 - severities)
    return pd.Series(data=W_new.values, index=W_aligned.index, name="W_sigmoid")


def build_graph(df):
    """Build a directed graph from the GENIE3 links DataFrame."""
    G = nx.DiGraph()
    for _, row in df.iterrows():
        G.add_edge(row["regulatoryGene"], row["targetGene"], weight=row["weight"])
    return G


def apply_perturbation_graph(G, gene, am_scores, k=10, t=0.5):
    """Apply a perturbation to the graph by modifying the weights of edges
    originating from the perturbed gene based on its AM scores."""
    H = G.copy()
    severity = logistic_transform(am_scores[gene], k, t)
    scale = 1 - severity
    for u, v, data in H.out_edges(gene, data=True):
        data["weight"] *= scale
    return H


def propagate_impacts_graph(G, initial_impacts, steps=2):
    """Propagate impacts through the graph for a given number of steps."""
    nodes = list(G.nodes)
    idx = {g: i for i, g in enumerate(nodes)}
    adj = nx.to_numpy_array(G, nodelist=nodes, weight="weight")
    vec = np.zeros(len(nodes))
    for g, val in initial_impacts.items():
        vec[idx[g]] = val
    for _ in range(steps):
        vec = vec @ adj
    return {nodes[i]: vec[i] for i in range(len(nodes))}


def main():
    # Build base graph
    G = build_graph(network_df)

    results = []
    for cell_id, group in am_df.groupby("sample"):
        # Build per-cell AM scores
        am_scores_cell = dict(zip(group["gene"], group["score"]))
        genes = list(am_scores_cell.keys())
        # Apply all perturbations sequentially
        H = G.copy()
        for g in genes:
            H = apply_perturbation_graph(H, g, am_scores_cell, k=k, t=t)
        # Build initial impact vector
        initial_impacts = {
            g: logistic_transform(am_scores_cell[g], k, t) for g in genes
        }
        # Propagate impacts
        impacts = propagate_impacts_graph(H, initial_impacts, steps=steps)
        # Record results
        for affected, val in impacts.items():
            results.append(
                {
                    "cell_id": cell_id,
                    "perturbed_genes": ",".join(genes),
                    "affected_gene": affected,
                    "impact": val,
                }
            )

    # Write output table
    out_df = pd.DataFrame(results)
    out_df.to_csv(scores_out_file, sep="\t", index=False)


if __name__ == "__main__":
    main()
