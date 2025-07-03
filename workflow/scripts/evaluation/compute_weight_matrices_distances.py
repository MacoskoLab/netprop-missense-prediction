"""
Evaluate all weight matrices against each other using multiple distance metrics.
Outputs a TSV file with pairwise comparisons.
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.script import Snakemake

    snakemake: Snakemake
    snakemake = None  # type: ignore


import os
import re
from itertools import combinations

import pandas as pd
from scipy.spatial.distance import euclidean
from scipy.stats import energy_distance, wasserstein_distance


def compute_distance_metrics(matrix1, matrix2, name1, name2):
    """
    Compute distance metrics between two weight matrices.

    Matrices must have exactly the same genes for both rows and columns,
    and contain no NaN values.

    Returns:
        dict: Dictionary with euclidean, wasserstein, and energy distances
    """
    matrix1 = matrix1.fillna(0)
    matrix2 = matrix2.fillna(0)

    print(f"Comparing matrices: {name1} vs {name2}")

    # Check that matrices have exactly the same row indices (genes)
    if not matrix1.index.equals(matrix2.index):
        raise ValueError(
            f"Matrices have different row genes. "
            f"{name1} has {len(matrix1.index)} genes, {name2} has {len(matrix2.index)} genes. "
            f"Row genes must be identical."
        )

    # Check that matrices have exactly the same column indices (genes)
    if not matrix1.columns.equals(matrix2.columns):
        raise ValueError(
            f"Matrices have different column genes. "
            f"{name1} has {len(matrix1.columns)} genes, {name2} has {len(matrix2.columns)} genes. "
            f"Column genes must be identical."
        )

    # Check for NaN values in matrix1
    if matrix1.isna().any().any():
        print(f"WARNING: {name1} contains NaN values.")

    # Check for NaN values in matrix2
    if matrix2.isna().any().any():
        print(f"WARNING: {name2} contains NaN values.")

    # Flatten matrices for distance computation
    m1_flat = matrix1.values.flatten()
    m2_flat = matrix2.values.flatten()

    # Compute distance metrics
    euclidean_dist = euclidean(m1_flat, m2_flat)
    wasserstein_dist = wasserstein_distance(m1_flat, m2_flat)
    energy_dist = energy_distance(m1_flat, m2_flat)

    return {
        "euclidean_distance": euclidean_dist,
        "wasserstein_distance": wasserstein_dist,
        "energy_distance": energy_dist,
    }


def load_matrix(filepath) -> pd.DataFrame:
    """
    Read a weight matrix from a file.

    Args:
        filepath (str): Path to the file containing the weight matrix.

    Returns:
        pd.DataFrame: DataFrame containing the weight matrix.
    """
    if filepath.endswith(".h5"):
        return pd.read_hdf(filepath, key="weight_matrix", mode="r")  # type: ignore
    elif filepath.endswith(".tsv") or filepath.endswith(".txt"):
        return pd.read_csv(filepath, sep="\t", index_col=0, header=0)
    else:
        raise ValueError(
            f"Unsupported file format for {filepath}. Only .h5, .tsv, or .txt are supported."
        )


def compare_matrices(matrices):
    """
    Compare a set of pairs of matrices and compute distance metrics.

    Args:
        matrices (dict): Dictionary with matrix names as keys and DataFrames as values.

    Returns:
        list: List of dictionaries with comparison results.
    """
    results = []
    for (name1, matrix1), (name2, matrix2) in combinations(matrices.items(), 2):
        metrics = compute_distance_metrics(matrix1, matrix2, name1, name2)
        results.append({"matrix1": name1, "matrix2": name2, **metrics})
    return pd.DataFrame(results)


def main():
    # Get input files from snakemake
    real_perturbed = load_matrix(snakemake.input["real_perturbed_matrix"])
    real_unperturbed = load_matrix(snakemake.input["real_unperturbed_matrix"])
    matrix_combinations = pd.read_csv(snakemake.input["matrix_combinations"], sep="\t")
    predicted_perturbed_matrices = []
    print(
        f"Looking for predicted matrices in: {snakemake.input['predicted_perturbed_matrices']}"
    )

    for matrix_fp in snakemake.input["predicted_perturbed_matrices"]:
        # Extract combination ID from filename using regex
        filename = os.path.basename(matrix_fp)
        print(f"Processing file: {filename}")
        match = re.search(r"predicted_perturbed_weights_(\d+)\.", filename)
        if match:
            combination_id = int(match.group(1))
            matrix = load_matrix(matrix_fp)
            print(f"Loaded matrix for combination {combination_id}")
            predicted_perturbed_matrices.append((matrix, combination_id))
        else:
            raise ValueError(
                f"Could not extract combination ID from filename: {filename}"
            )

    print(f"Found {len(predicted_perturbed_matrices)} predicted matrices")
    if len(predicted_perturbed_matrices) == 0:
        print("ERROR: No predicted matrices found!")
        return

    output_fp = snakemake.output.results

    comparison_results = []

    for predicted_matrix, combination_id in predicted_perturbed_matrices:
        print(f"Processing combination ID: {combination_id}")
        # Load matrices
        matrices = {
            "real_perturbed": real_perturbed,
            "real_unperturbed": real_unperturbed,
            "predicted_perturbed": predicted_matrix,
        }

        # Compare matrices
        comparison = compare_matrices(matrices)
        comparison["combination_id"] = combination_id
        comparison_results.append(comparison)

    results_df = pd.concat(comparison_results, ignore_index=True)

    # Add results to grid params
    results_df = results_df.merge(matrix_combinations, on="combination_id", how="left")

    results_df.to_csv(output_fp, sep="\t", index=False)
    print(f"Results saved to {output_fp}")
    print("All pairwise comparisons completed.")


if __name__ == "__main__":
    main()
