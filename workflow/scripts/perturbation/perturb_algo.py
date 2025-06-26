from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.script import Snakemake

    snakemake: Snakemake
    snakemake = None  # type: ignore

import numpy as np
import pandas as pd


def threshold_transform_matrix(
    weight_matrix: pd.DataFrame, score_series: pd.Series, threshold: float
) -> pd.DataFrame:
    """
    For each gene with score >= threshold, scale outgoing edge weights (rows).
    """
    weight_matrix_out = weight_matrix.copy()

    # Create a series aligned with the weight matrix index, filling missing genes with 0
    aligned_scores = score_series.reindex(weight_matrix.index, fill_value=0.0)

    # Vectorized computation: handle NaN values and compute fractions
    valid_mask = ~np.isnan(aligned_scores) & (aligned_scores >= threshold)
    fractions = np.zeros_like(aligned_scores, dtype=float)
    fractions[valid_mask] = np.clip(
        (aligned_scores[valid_mask] - threshold) / (1.0 - threshold), 0.0, 1.0
    )

    # Vectorized scaling: multiply each row by (1 - fraction)
    scaling_factors = 1.0 - fractions
    weight_matrix_out = weight_matrix_out.multiply(scaling_factors, axis="index")

    print("threshold_transform_matrix finished", flush=True)
    return weight_matrix_out


def sigmoid_transform_matrix(
    weight_matrix: pd.DataFrame,
    score_series: pd.Series,
    steepness: float,
    midpoint: float,
) -> pd.DataFrame:
    """
    For each gene, apply sigmoid-based scaling to outgoing edge weights (rows).
    """
    weight_matrix_out = weight_matrix.copy()

    # Create a series aligned with the weight matrix index, filling missing genes with 0
    aligned_scores = score_series.reindex(weight_matrix.index, fill_value=0.0)

    # Vectorized sigmoid computation: handle NaN values
    valid_mask = ~np.isnan(aligned_scores)
    fractions = np.zeros_like(aligned_scores, dtype=float)
    fractions[valid_mask] = 1.0 / (
        1.0 + np.exp(-steepness * (aligned_scores[valid_mask] - midpoint))
    )

    # Vectorized scaling: multiply each row by (1 - fraction)
    scaling_factors = 1.0 - fractions
    weight_matrix_out = weight_matrix_out.multiply(scaling_factors, axis="index")

    print("sigmoid_transform_matrix finished", flush=True)
    return weight_matrix_out


def propagate_effects_matrix(
    weight_matrix: pd.DataFrame, score_series: pd.Series, steps: int
) -> pd.Series:
    """
    Propagate initial scores through matrix multiplication.
    Returns a series of cumulative perturbation per gene.
    """
    # Ensure all genes in the matrix are represented in the propagated series
    all_genes = weight_matrix.index.union(weight_matrix.columns)
    propagated = pd.Series({gene: score_series.get(gene, 0.0) for gene in all_genes})

    for _ in range(steps):
        # Matrix multiplication: propagated values flow through edges
        # weight_matrix.T because we want to multiply by the transpose
        # (target genes as rows, source genes as columns)
        propagated_next = weight_matrix.T.dot(
            propagated.reindex(weight_matrix.index, fill_value=0.0)
        )
        propagated = propagated.add(propagated_next, fill_value=0.0)

    print("propagate_effects_matrix finished", flush=True)
    return propagated


def update_matrix_weights_with_propagation(
    weight_matrix_original: pd.DataFrame, propagated: pd.Series
) -> pd.DataFrame:
    """
    Scale each gene's outgoing edges (rows) by (1 - propagated[gene]).
    """
    weight_matrix_out = weight_matrix_original.copy()

    # Create a series aligned with the weight matrix index, filling missing genes with 0
    aligned_propagated = propagated.reindex(
        weight_matrix_original.index, fill_value=0.0
    )

    # Vectorized scaling: multiply each row by (1 - propagated_value)
    scaling_factors = 1.0 - aligned_propagated
    weight_matrix_out = weight_matrix_out.multiply(scaling_factors, axis="index")

    print("update_matrix_weights_with_propagation finished", flush=True)
    return weight_matrix_out


def main():
    # Access parameters directly from config
    perturb_config = snakemake.config["perturbation_algorithm"]

    steps = perturb_config.get("steps")
    method = perturb_config.get("score_transform")

    # Load file paths
    genie3_weights_file_path = snakemake.input["genie3_weights"]
    am_scores_path = snakemake.input["am_scores"]
    perturb_list_path = snakemake.input["perturbations_list"]
    output_path = snakemake.output["perturbed_weights"]

    # Load weight matrix from TSV
    weight_matrix = pd.read_csv(genie3_weights_file_path, sep="\t", index_col=0)

    am_df = pd.read_csv(am_scores_path, sep="\t")
    perturb_df = pd.read_csv(perturb_list_path, sep="\t")

    # Filter AM scores to only those gene-variant pairs
    sel = am_df.merge(perturb_df, on=["gene", "variant"])
    # Build score series indexed by gene; exact variants, no aggregation
    S_series = sel.set_index("gene")["score"]

    print("Beginning perturbation algorithm...", flush=True)

    if method == "threshold":
        threshold_value = perturb_config.get("threshold")
        weight_matrix_scaled = threshold_transform_matrix(
            weight_matrix, S_series, threshold_value
        )
    elif method == "sigmoid":
        steepness = perturb_config.get("steepness")
        midpoint = perturb_config.get("midpoint")
        weight_matrix_scaled = sigmoid_transform_matrix(
            weight_matrix, S_series, steepness, midpoint
        )
    else:
        raise ValueError(f"Unknown method: {method}")

    propagated_series = propagate_effects_matrix(weight_matrix_scaled, S_series, steps)

    weight_matrix_final = update_matrix_weights_with_propagation(
        weight_matrix, propagated_series
    )

    # Export weight matrix
    weight_matrix_final.to_hdf(output_path, key="weight_matrix", mode="w")
    print(f"Weight matrix written to {output_path}", flush=True)


if __name__ == "__main__":
    main()
