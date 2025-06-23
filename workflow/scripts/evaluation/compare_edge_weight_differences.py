from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.script import Snakemake

    snakemake: Snakemake
    snakemake = None  # type: ignore

import numpy as np
import pandas as pd
from scipy.spatial.distance import cosine
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error


def compute_edge_weight_differences(
    control_df: pd.DataFrame, perturbed_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Compute edge weight differences: Δ = w_perturbed - w_control

    Args:
        control_df: Control network with columns: regulatoryGene, targetGene, weight
        perturbed_df: Perturbed network with columns: regulatoryGene, targetGene, weight

    Returns:
        DataFrame with edge weight differences
    """
    # Merge on regulatory and target genes
    merged = pd.merge(
        control_df,
        perturbed_df,
        on=["regulatoryGene", "targetGene"],
        suffixes=("_control", "_perturbed"),
        how="inner",  # Only include edges present in both networks
    )

    # Compute differences
    merged["weight_diff"] = merged["weight_perturbed"] - merged["weight_control"]

    # Select relevant columns
    result = merged[
        [
            "regulatoryGene",
            "targetGene",
            "weight_control",
            "weight_perturbed",
            "weight_diff",
        ]
    ].copy()

    return result


def compare_weight_differences(
    real_diff_df: pd.DataFrame, pred_diff_df: pd.DataFrame
) -> dict:
    """
    Compare real vs predicted edge weight differences.

    Args:
        real_diff_df: Real edge weight differences
        pred_diff_df: Predicted edge weight differences

    Returns:
        Dictionary with comparison metrics
    """
    # Merge on regulatory and target genes to align edges
    merged = pd.merge(
        real_diff_df[["regulatoryGene", "targetGene", "weight_diff"]],
        pred_diff_df[["regulatoryGene", "targetGene", "weight_diff"]],
        on=["regulatoryGene", "targetGene"],
        suffixes=("_real", "_pred"),
        how="inner",
    )

    print(f"Comparing {len(merged)} common edges")

    if len(merged) == 0:
        print("Warning: No common edges found between real and predicted differences")
        return {
            "n_edges": 0,
            "pearson_r": np.nan,
            "pearson_p": np.nan,
            "spearman_r": np.nan,
            "spearman_p": np.nan,
            "mse": np.nan,
            "mae": np.nan,
            "cosine_distance": np.nan,
        }

    real_diffs = merged["weight_diff_real"].values
    pred_diffs = merged["weight_diff_pred"].values

    # Remove any NaN values
    valid_mask = ~(np.isnan(real_diffs) | np.isnan(pred_diffs))  # type: ignore
    real_diffs = real_diffs[valid_mask]
    pred_diffs = pred_diffs[valid_mask]

    if len(real_diffs) == 0:
        print("Warning: No valid (non-NaN) edge differences found")
        return {
            "n_edges": 0,
            "pearson_r": np.nan,
            "pearson_p": np.nan,
            "spearman_r": np.nan,
            "spearman_p": np.nan,
            "mse": np.nan,
            "mae": np.nan,
            "cosine_distance": np.nan,
        }

    # Compute correlation metrics
    try:
        pearson_r, pearson_p = pearsonr(real_diffs, pred_diffs)
    except:
        pearson_r, pearson_p = np.nan, np.nan

    try:
        spearman_r, spearman_p = spearmanr(real_diffs, pred_diffs)
    except:
        spearman_r, spearman_p = np.nan, np.nan

    # Compute error metrics
    mse = mean_squared_error(real_diffs, pred_diffs)
    mae = mean_absolute_error(real_diffs, pred_diffs)

    # Compute cosine distance
    try:
        cosine_dist = cosine(real_diffs, pred_diffs)
    except:
        cosine_dist = np.nan

    return {
        "n_edges": len(real_diffs),
        "pearson_r": pearson_r,
        "pearson_p": pearson_p,
        "spearman_r": spearman_r,
        "spearman_p": spearman_p,
        "mse": mse,
        "mae": mae,
        "cosine_distance": cosine_dist,
    }


def main():
    """Main function to compare edge weight differences."""

    # Load input data
    print("Loading edge weight data...")
    real_control_weights = pd.read_csv(snakemake.input.real_control_weights, sep="\t")
    real_perturbed_weights = pd.read_csv(
        snakemake.input.real_perturbed_weights, sep="\t"
    )
    pred_control_weights = pd.read_csv(
        snakemake.input.pred_control_weights, sep="\t"
    )  # This is the original control network
    pred_perturbed_weights = pd.read_csv(
        snakemake.input.pred_perturbed_weights, sep="\t"
    )

    print(f"Real control weights: {len(real_control_weights)} edges")
    print(f"Real perturbed weights: {len(real_perturbed_weights)} edges")
    print(f"Predicted control weights: {len(pred_control_weights)} edges")
    print(f"Predicted perturbed weights: {len(pred_perturbed_weights)} edges")

    # Compute edge weight differences
    print("Computing real edge weight differences...")
    real_diff_df = compute_edge_weight_differences(
        real_control_weights, real_perturbed_weights
    )

    print("Computing predicted edge weight differences...")
    pred_diff_df = compute_edge_weight_differences(
        pred_control_weights, pred_perturbed_weights
    )

    print(f"Real differences: {len(real_diff_df)} edges")
    print(f"Predicted differences: {len(pred_diff_df)} edges")

    # Compare the differences
    print("Comparing edge weight differences...")
    comparison_metrics = compare_weight_differences(real_diff_df, pred_diff_df)

    # Save edge weight differences
    real_diff_df.to_csv(snakemake.output.real_edge_diffs, sep="\t", index=False)
    pred_diff_df.to_csv(snakemake.output.pred_edge_diffs, sep="\t", index=False)

    # Save comparison results
    comparison_df = pd.DataFrame([comparison_metrics])
    comparison_df.to_csv(snakemake.output.comparison_results, sep="\t", index=False)

    print("Edge weight difference comparison completed")
    print(f"Real edge differences saved to: {snakemake.output.real_edge_diffs}")
    print(f"Predicted edge differences saved to: {snakemake.output.pred_edge_diffs}")
    print(f"Comparison results saved to: {snakemake.output.comparison_results}")

    # Print key metrics
    print("\n=== COMPARISON RESULTS ===")
    print(f"Number of common edges: {comparison_metrics['n_edges']}")
    print(
        f"Pearson correlation: {comparison_metrics['pearson_r']:.4f} (p={comparison_metrics['pearson_p']:.4e})"
    )
    print(
        f"Spearman correlation: {comparison_metrics['spearman_r']:.4f} (p={comparison_metrics['spearman_p']:.4e})"
    )
    print(f"Mean Squared Error: {comparison_metrics['mse']:.4f}")
    print(f"Mean Absolute Error: {comparison_metrics['mae']:.4f}")
    print(f"Cosine distance: {comparison_metrics['cosine_distance']:.4f}")


if __name__ == "__main__":
    main()
