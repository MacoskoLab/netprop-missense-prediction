from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.script import Snakemake

    snakemake: Snakemake
    snakemake = None  # type: ignore

import warnings

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")


def load_expression_data(expr_path: str) -> pd.DataFrame:
    """Load expression matrix from file."""
    return pd.read_csv(expr_path, sep="\t", index_col=0)


def compute_edge_weights_from_expression(
    expr_df: pd.DataFrame, network_df: pd.DataFrame
) -> pd.DataFrame:
    """
    For each edge in the network, compute edge weight by fitting linear regression:
    target_expr ~ regulator_expr

    Args:
        expr_df: Expression matrix with genes as rows and samples as columns
        network_df: Network dataframe with columns: regulatoryGene, targetGene, weight

    Returns:
        DataFrame with columns: regulatoryGene, targetGene, weight
    """
    # blend parameter for original and regression weights
    alpha = snakemake.params.get("blend_alpha", 0.5)
    results = []

    # Convert expression matrix so genes are columns (features) and samples are rows
    expr_transposed = expr_df.T

    for _, row in network_df.iterrows():
        regulator = row["regulatoryGene"]
        target = row["targetGene"]

        # Check if both genes exist in expression data
        if (
            regulator not in expr_transposed.columns
            or target not in expr_transposed.columns
        ):
            # Use original weight if genes not found
            results.append(
                {
                    "regulatoryGene": regulator,
                    "targetGene": target,
                    "weight": row["weight"],
                }
            )
            continue

        # Get expression vectors
        X = expr_transposed[regulator].values.reshape(-1, 1)  # regulator expression
        y = expr_transposed[target].values  # target expression

        # Fit linear regression
        try:
            reg = LinearRegression()
            reg.fit(X, y)

            # Use coefficient and blend with original weight
            coef = reg.coef_[0]
            original = row["weight"]
            weight = alpha * original + (1.0 - alpha) * coef

            results.append(
                {"regulatoryGene": regulator, "targetGene": target, "weight": weight}
            )

        except Exception as e:
            print(f"Error fitting regression for {regulator} -> {target}: {e}")
            # Use original weight as fallback
            results.append(
                {
                    "regulatoryGene": regulator,
                    "targetGene": target,
                    "weight": row["weight"],
                }
            )

    return pd.DataFrame(results)


def main():
    """Main function to compute real perturbed edge weights using control topology."""

    # Load inputs
    control_network = pd.read_csv(snakemake.input.control_network, sep="\t")
    perturbed_expr = load_expression_data(snakemake.input.perturbed_expr)
    control_expr = load_expression_data(snakemake.input.control_expr)

    print(f"Control network has {len(control_network)} edges")
    print(
        f"Perturbed expression data: {perturbed_expr.shape[0]} genes, {perturbed_expr.shape[1]} samples"
    )
    print(
        f"Control expression data: {control_expr.shape[0]} genes, {control_expr.shape[1]} samples"
    )

    # Keep the original GENIE3 control weights unchanged
    print("Using original GENIE3 control weights (no recomputation)...")
    control_weights_df = control_network.copy()

    # Only compute new edge weights for perturbed data using control topology
    print("Computing perturbed edge weights using control topology...")
    perturbed_weights_df = compute_edge_weights_from_expression(
        perturbed_expr, control_network
    )

    # Save results
    perturbed_weights_df.to_csv(
        snakemake.output.perturbed_weights, sep="\t", index=False
    )
    control_weights_df.to_csv(snakemake.output.control_weights, sep="\t", index=False)

    print("Edge weight computation completed")
    print(
        f"Perturbed weights (linear regression) saved to: {snakemake.output.perturbed_weights}"
    )
    print(
        f"Control weights (original GENIE3) saved to: {snakemake.output.control_weights}"
    )
    print("Now comparing: LinearRegression_perturbed vs GENIE3_control")


if __name__ == "__main__":
    main()
