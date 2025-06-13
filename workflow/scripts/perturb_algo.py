from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.script import Snakemake

    snakemake: Snakemake
    snakemake = None  # type: ignore

import pandas as pd
import numpy as np
import networkx as nx


def build_graph(network_df: pd.DataFrame) -> nx.DiGraph:
    """
    Construct a directed graph with 'weight' attribute on edges.
    """
    G = nx.DiGraph()
    for source_gene, target_gene, weight in zip(
        network_df["regulatoryGene"], network_df["targetGene"], network_df["weight"]
    ):
        G.add_edge(source_gene, target_gene, weight=weight)
    return G


def threshold_transform_graph(
    G: nx.DiGraph, score_series: pd.Series, threshold: float
) -> nx.DiGraph:
    """
    For each node with score >= threshold, scale outgoing edge weights.
    """
    G_out = G.copy()
    for node, score in score_series.items():
        if np.isnan(score) or score < threshold:
            continue
        frac = np.clip((score - threshold) / (1.0 - threshold), 0.0, 1.0)
        for _, tgt, data in G_out.out_edges(node, data=True):  # type: ignore
            data["weight"] *= 1.0 - frac
    return G_out  # type: ignore


def sigmoid_transform_graph(
    G: nx.DiGraph, score_series: pd.Series, steepness: float, midpoint: float
) -> nx.DiGraph:
    """
    For each node, apply sigmoid-based scaling to outgoing edge weights.
    """
    G_out = G.copy()
    for node, score in score_series.items():
        if np.isnan(score):
            continue
        frac = 1.0 / (1.0 + np.exp(-steepness * (score - midpoint)))
        for _, tgt, data in G_out.out_edges(node, data=True):  # type: ignore
            data["weight"] *= 1.0 - frac
    return G_out  # type: ignore


def propagate_effects_graph(
    G: nx.DiGraph, score_series: pd.Series, steps: int
) -> pd.Series:
    """
    Propagate initial scores through graph adjacency weights.
    Returns a series of cumulative perturbation per node.
    """
    nodes = list(G.nodes)
    propagated = pd.Series({n: score_series.get(n, 0.0) for n in nodes})
    for _ in range(steps):
        propagated_next = pd.Series(0.0, index=nodes)
        for source_gene, target_gene, data in G.edges(data=True):
            propagated_next[target_gene] += propagated[source_gene] * data["weight"]
        propagated += propagated_next
    return propagated


def update_graph_weights_with_propagation(
    G_original: nx.DiGraph, propagated: pd.Series
) -> nx.DiGraph:
    """
    Scale each node's outgoing edges by (1 - propagated[node]).
    """
    G_out = G_original.copy()
    for node in G_out.nodes:
        # if not in path we are just multiplying by 1
        factor = 1.0 - propagated.get(node, 0.0)
        for _, tgt, data in G_out.out_edges(node, data=True):  # type: ignore
            data["weight"] *= factor
    return G_out  # type: ignore


def save_graph_results(G: nx.DiGraph, output_path: str):
    """
    Write edge list with weights to file.
    """
    rows = [(src, tgt, data["weight"]) for src, tgt, data in G.edges(data=True)]
    out_df = pd.DataFrame(rows, columns=["regulatoryGene", "targetGene", "weight"])
    out_df.to_csv(output_path, sep="\t", index=False, header=True)


def main():
    # Load parameters
    steps = snakemake.params.get("steps", 2)
    steepness = snakemake.params.get("steepness", 10.0)
    midpoint = snakemake.params.get("midpoint", 0.5)
    threshold_value = snakemake.params.get("threshold", 0.6)
    method = (
        snakemake.params.get("pathogenicity_score_transform_method", "threshold")
        .strip()
        .lower()
    )

    # Load file paths
    genie3_links_path = snakemake.input["genie3_links"]
    am_scores_path = snakemake.input["am_scores"]
    perturb_list_path = snakemake.input["perturbations_list"]
    output_path = snakemake.output["scores"]

    # Load data
    network_df = pd.read_csv(genie3_links_path, sep="\t")
    am_df = pd.read_csv(am_scores_path, sep="\t")
    perturb_df = pd.read_csv(perturb_list_path, sep="\t")
    # Filter AM scores to only those gene-variant pairs
    sel = am_df.merge(perturb_df, on=["gene", "variant"])
    # Build score series indexed by gene; exact variants, no aggregation
    S_series = sel.set_index("gene")["score"]

    G = build_graph(network_df)
    if method == "threshold":
        G_scaled = threshold_transform_graph(G, S_series, threshold_value)
    elif method == "sigmoid":
        G_scaled = sigmoid_transform_graph(G, S_series, steepness, midpoint)
    else:
        raise ValueError(f"Unknown method: {method}")

    propagated_series = propagate_effects_graph(G_scaled, S_series, steps)
    G_final = update_graph_weights_with_propagation(G, propagated_series)
    save_graph_results(G_final, output_path)


main()
