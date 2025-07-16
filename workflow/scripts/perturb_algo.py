from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.script import Snakemake

    snakemake: Snakemake
    snakemake = None  # type: ignore

import numpy as np
import pandas as pd


def threshold_transform_matrix(
    weight_matrix: pd.DataFrame, gene: str, score: float, threshold: float
) -> pd.DataFrame:
    """
    For the specified gene with score >= threshold, scale its outgoing edge weights (row).
    """
    weight_matrix_out = weight_matrix.copy()

    # Only transform if score meets threshold
    if score >= threshold:
        fraction = np.clip((score - threshold) / (1.0 - threshold), 0.0, 1.0)
        scaling_factor = 1.0 - fraction
        weight_matrix_out.loc[gene] *= scaling_factor

    return weight_matrix_out


def sigmoid_transform_matrix(
    weight_matrix: pd.DataFrame,
    gene: str,
    score: float,
    steepness: float,
    midpoint: float,
) -> pd.DataFrame:
    """
    For the specified gene, apply sigmoid-based scaling to its outgoing edge weights (row).
    """
    weight_matrix_out = weight_matrix.copy()

    fraction = 1.0 / (1.0 + np.exp(-steepness * (score - midpoint)))
    scaling_factor = 1.0 - fraction
    weight_matrix_out.loc[gene] *= scaling_factor

    return weight_matrix_out


def propagate_effects_matrix(
    weight_matrix: pd.DataFrame, gene: str, score: float, steps: int
) -> pd.Series:
    """
    Propagate initial score from a single gene through matrix multiplication.
    Returns a series of cumulative perturbation per gene.
    """
    # Initialize propagated values for all genes in the weight matrix
    propagated = pd.Series(0.0, index=weight_matrix.index)

    # Set initial score for the perturbed gene
    propagated[gene] = score

    if debug:
        print(f"Initial propagation - {gene}: {score}", flush=True)
        print(f"Step 0: Total propagated effect = {propagated.sum():.6f}", flush=True)
        print(f"Step 0: Non-zero genes = {(propagated != 0).sum()}", flush=True)

    for step in range(steps):
        propagated_next = weight_matrix.T.dot(propagated)
        propagated = propagated.add(propagated_next, fill_value=0.0)

        if debug:
            print(
                f"Step {step + 1}: Total propagated effect = {propagated.sum():.6f}",
                flush=True,
            )
            print(
                f"Step {step + 1}: Non-zero genes = {(propagated != 0).sum()}",
                flush=True,
            )
            print(
                f"Step {step + 1}: Max effect = {propagated.max():.6f}, Min effect = {propagated.min():.6f}",
                flush=True,
            )

    if debug:
        print(
            f"Final propagation complete. Affected genes: {(propagated != 0).sum()}/{len(propagated)}",
            flush=True,
        )

    return propagated


def update_matrix_weights_with_propagation(
    weight_matrix_original: pd.DataFrame, propagated: pd.Series
) -> pd.DataFrame:
    """
    Scale each gene's outgoing edges (rows) by (1 - propagated[gene]).
    """
    if debug:
        print(
            f"Updating weight matrix with propagated effects. Affected genes: {(propagated != 0).sum()}/{len(propagated)}",
            flush=True,
        )

    weight_matrix_out = weight_matrix_original.copy()
    scaling_factors = 1.0 - propagated
    weight_matrix_out = weight_matrix_out.multiply(scaling_factors, axis="index")

    if debug:
        print(
            f"Weight matrix updated. First 10 edge weights for {propagated.index[0]}: {weight_matrix_out.loc[propagated.index[0]].values[:10]}",
            flush=True,
        )

    return weight_matrix_out


def main():
    # Handle debug flag
    global debug
    debug = snakemake.config["perturbation_algorithm"].get("debug", False)
    print(f"Debug mode: {debug}", flush=True)

    # Handle combination ID
    combination_id = int(snakemake.params["combination_id"])
    combinations_df = pd.read_csv(snakemake.input["combinations"], sep="\t")
    combination_row = combinations_df[
        combinations_df["combination_id"] == combination_id
    ].iloc[0]

    # Load files
    print("Loading weight matrix...", flush=True)
    weight_matrix = pd.read_csv(
        snakemake.input["genie3_weights"], sep="\t", index_col=0
    )
    print(f"Weight matrix shape: {weight_matrix.shape}", flush=True)
    am_df = pd.read_csv(snakemake.input["am_scores"], sep="\t")
    perturb_df = pd.read_csv(snakemake.input["perturbations_list"], sep="\t")

    perturbed_gene = perturb_df["gene"].iloc[0]  # Expect exactly one gene-variant pair
    perturbed_score = am_df.loc[am_df["gene"] == perturbed_gene, "score"].iloc[0]

    if debug:
        print(
            f"Perturbing gene: {perturbed_gene} with AlphaMissense score: {perturbed_score}"
        )

    # Sanity check!
    if perturbed_gene not in weight_matrix.index:
        raise ValueError(
            f"Perturbed gene '{perturbed_gene}' not found in weight matrix index."
        )

    if debug:
        print(
            f"First 10 edge weights for {perturbed_gene} before scaling: {weight_matrix.loc[perturbed_gene].values[:10]}"
        )

    if combination_row["score_transform"] == "threshold":
        weight_matrix_scaled = threshold_transform_matrix(
            weight_matrix,
            perturbed_gene,
            perturbed_score,
            combination_row["threshold"],
        )
    elif combination_row["score_transform"] == "sigmoid":
        weight_matrix_scaled = sigmoid_transform_matrix(
            weight_matrix,
            perturbed_gene,
            perturbed_score,
            combination_row["steepness"],
            combination_row["midpoint"],
        )
    else:
        raise ValueError(f"Unknown method: {combination_row['score_transform']}")

    if debug:
        print(
            f"First 10 edge weights for {perturbed_gene} after scaling: {weight_matrix_scaled.loc[perturbed_gene].values[:10]}"
        )

    propagated_series = propagate_effects_matrix(
        weight_matrix_scaled,
        perturbed_gene,
        perturbed_score,
        int(combination_row["steps"]),
    )

    weight_matrix_final = update_matrix_weights_with_propagation(
        weight_matrix, propagated_series
    )

    # Export weight matrix
    weight_matrix_final.to_hdf(
        snakemake.output["perturbed_weights"], key="weight_matrix", mode="w"
    )


if __name__ == "__main__":
    main()
