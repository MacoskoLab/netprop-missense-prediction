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
    For the specified gene with score >= threshold, scale the entire weight matrix.
    """
    weight_matrix_out = weight_matrix.copy()

    # Only transform if score meets threshold
    if score >= threshold:
        fraction = np.clip((score - threshold) / (1.0 - threshold), 0.0, 1.0)
        scaling_factor = 1.0 - fraction
        weight_matrix_out *= scaling_factor

    return weight_matrix_out


def sigmoid_transform_matrix(
    weight_matrix: pd.DataFrame,
    gene: str,
    score: float,
    steepness: float,
    midpoint: float,
) -> pd.DataFrame:
    """
    Apply sigmoid scaling to the entire weight matrix based on the gene score.
    """
    weight_matrix_out = weight_matrix.copy()

    fraction = 1.0 / (1.0 + np.exp(-steepness * (score - midpoint)))
    scaling_factor = 1.0 - fraction
    weight_matrix_out *= scaling_factor

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
        print(
            f"Weight matrix indices check - Rows: {len(weight_matrix.index)}, Cols: {len(weight_matrix.columns)}",
            flush=True,
        )
        print(
            f"Row and column indices identical: {weight_matrix.index.equals(weight_matrix.columns)}",
            flush=True,
        )
        print(f"Initial propagation - {gene}: {score}", flush=True)
        print(f"Step 0: Total propagated effect = {propagated.sum():.6f}", flush=True)
        print(f"Step 0: Non-zero genes = {(propagated != 0).sum()}", flush=True)

    for step in range(steps):
        propagated_next = weight_matrix.T.dot(propagated)
        # Ensure propagated_next has the same index as propagated to avoid index expansion
        propagated_next = propagated_next.reindex(propagated.index, fill_value=0.0)
        propagated = propagated + propagated_next

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


def compare_matrices(
    predicted_matrix: pd.DataFrame,
    real_matrix: pd.DataFrame,
    propagated_series: pd.Series,
) -> None:
    """
    Compare predicted and real perturbed matrices, separating analysis for
    genes affected vs unaffected by the propagation algorithm.
    """
    print("\n" + "=" * 60, flush=True)
    print("MATRIX COMPARISON ANALYSIS", flush=True)
    print("=" * 60, flush=True)

    # Ensure matrices have the same shape and indices
    common_genes = predicted_matrix.index.intersection(real_matrix.index)
    common_targets = predicted_matrix.columns.intersection(real_matrix.columns)

    if len(common_genes) == 0 or len(common_targets) == 0:
        print(
            "ERROR: No common genes/targets between predicted and real matrices!",
            flush=True,
        )
        return

    # Subset both matrices to common genes and targets
    pred_subset = predicted_matrix.loc[common_genes, common_targets]
    real_subset = real_matrix.loc[common_genes, common_targets]

    print(
        f"Comparing matrices with {len(common_genes)} genes and {len(common_targets)} targets",
        flush=True,
    )

    # Identify affected genes (those with non-zero propagation effects)
    affected_genes = propagated_series[propagated_series != 0].index
    affected_genes = affected_genes.intersection(common_genes)
    unaffected_genes = common_genes.difference(affected_genes)

    print(f"Affected genes by propagation: {len(affected_genes)}", flush=True)
    print(f"Unaffected genes by propagation: {len(unaffected_genes)}", flush=True)

    # Calculate differences
    diff_matrix = pred_subset - real_subset

    # Analysis for AFFECTED genes
    if len(affected_genes) > 0:
        print(
            f"\n--- AFFECTED GENES ANALYSIS (n={len(affected_genes)}) ---", flush=True
        )
        affected_pred = pred_subset.loc[affected_genes]
        affected_real = real_subset.loc[affected_genes]
        affected_diff = diff_matrix.loc[affected_genes]

        # Summary statistics
        print(
            f"Predicted values - Mean: {affected_pred.values.mean():.6f}, Std: {affected_pred.values.std():.6f}",
            flush=True,
        )
        print(
            f"Real values - Mean: {affected_real.values.mean():.6f}, Std: {affected_real.values.std():.6f}",
            flush=True,
        )
        print(
            f"Differences - Mean: {affected_diff.values.mean():.6f}, Std: {affected_diff.values.std():.6f}",
            flush=True,
        )
        print(
            f"Differences - Min: {affected_diff.values.min():.6f}, Max: {affected_diff.values.max():.6f}",
            flush=True,
        )

        # Correlation
        correlation = np.corrcoef(
            affected_pred.values.flatten(), affected_real.values.flatten()
        )[0, 1]
        print(f"Correlation (predicted vs real): {correlation:.6f}", flush=True)

        # MAE and RMSE
        mae = np.mean(np.abs(affected_diff.values))
        rmse = np.sqrt(np.mean(affected_diff.values**2))
        print(f"Mean Absolute Error: {mae:.6f}", flush=True)
        print(f"Root Mean Square Error: {rmse:.6f}", flush=True)

        # Show most affected genes
        gene_level_effects = (
            propagated_series[affected_genes].abs().sort_values(ascending=False)
        )
        print(f"Top 5 most affected genes by propagation:", flush=True)
        for i, (gene, effect) in enumerate(gene_level_effects.head().items()):
            gene_mae = np.mean(np.abs(affected_diff.loc[gene].values))
            print(
                f"  {i+1}. {gene}: propagation_effect={effect:.6f}, MAE={gene_mae:.6f}",
                flush=True,
            )

    # Analysis for UNAFFECTED genes
    if len(unaffected_genes) > 0:
        print(
            f"\n--- UNAFFECTED GENES ANALYSIS (n={len(unaffected_genes)}) ---",
            flush=True,
        )
        unaffected_pred = pred_subset.loc[unaffected_genes]
        unaffected_real = real_subset.loc[unaffected_genes]
        unaffected_diff = diff_matrix.loc[unaffected_genes]

        # Summary statistics
        print(
            f"Predicted values - Mean: {unaffected_pred.values.mean():.6f}, Std: {unaffected_pred.values.std():.6f}",
            flush=True,
        )
        print(
            f"Real values - Mean: {unaffected_real.values.mean():.6f}, Std: {unaffected_real.values.std():.6f}",
            flush=True,
        )
        print(
            f"Differences - Mean: {unaffected_diff.values.mean():.6f}, Std: {unaffected_diff.values.std():.6f}",
            flush=True,
        )
        print(
            f"Differences - Min: {unaffected_diff.values.min():.6f}, Max: {unaffected_diff.values.max():.6f}",
            flush=True,
        )

        # Correlation
        correlation = np.corrcoef(
            unaffected_pred.values.flatten(), unaffected_real.values.flatten()
        )[0, 1]
        print(f"Correlation (predicted vs real): {correlation:.6f}", flush=True)

        # MAE and RMSE
        mae = np.mean(np.abs(unaffected_diff.values))
        rmse = np.sqrt(np.mean(unaffected_diff.values**2))
        print(f"Mean Absolute Error: {mae:.6f}", flush=True)
        print(f"Root Mean Square Error: {rmse:.6f}", flush=True)

        # Check if unaffected genes should theoretically be identical
        if mae > 1e-10:  # Small tolerance for floating point errors
            print(
                f"WARNING: Unaffected genes show differences > 1e-10. This may indicate an issue.",
                flush=True,
            )
        else:
            print(
                f"GOOD: Unaffected genes show minimal differences (< 1e-10), as expected.",
                flush=True,
            )

    print("=" * 60 + "\n", flush=True)


def main():
    # Handle debug flag
    global debug
    debug = snakemake.config["perturbation_algorithm"].get("debug", False)
    print(f"Debug mode: {debug}", flush=True)

    global compare_with_real
    compare_with_real = snakemake.config["perturbation_algorithm"].get(
        "compare_with_real", False
    )
    print(f"Compare with real perturbed matrix: {compare_with_real}", flush=True)
    if compare_with_real:
        real_perturbed_weights = pd.read_csv(
            snakemake.input["real_perturbed_weights"], sep="\t", index_col=0
        )

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
        assert weight_matrix_scaled.index.equals(
            weight_matrix.index
        ) and weight_matrix_scaled.columns.equals(
            weight_matrix.columns
        ), "INDICES MISMATCH AFTER SCALING"

    propagated_series = propagate_effects_matrix(
        weight_matrix_scaled,
        perturbed_gene,
        perturbed_score,
        int(combination_row["steps"]),
    )

    weight_matrix_final = update_matrix_weights_with_propagation(
        weight_matrix, propagated_series
    )

    assert weight_matrix_final.index.equals(
        weight_matrix.index
    ) and weight_matrix_final.columns.equals(
        weight_matrix.columns
    ), "INDICES MISMATCH AFTER PROPAGATION"

    # Compare with real perturbed matrix if available
    if compare_with_real:
        compare_matrices(weight_matrix_final, real_perturbed_weights, propagated_series)

    # Export weight matrix
    weight_matrix_final.to_hdf(
        snakemake.output["perturbed_weights"], key="weight_matrix", mode="w"
    )


if __name__ == "__main__":
    main()
