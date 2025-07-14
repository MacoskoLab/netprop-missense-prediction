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
from itertools import combinations

import h5py
import pandas as pd
from scipy.spatial.distance import euclidean
from scipy.stats import energy_distance, wasserstein_distance


def compute_distance_metrics(matrix1_fp, matrix2_fp, name1, name2):
    """
    Compute distance metrics between two weight matrices.

    Matrices must have exactly the same genes for both rows and columns,
    and contain no NaN values.

    Returns:
        dict: Dictionary with euclidean, wasserstein, and energy distances
    """
    matrix1 = load_matrix(matrix1_fp, name1).fillna(0)
    matrix2 = load_matrix(matrix2_fp, name2).fillna(0)

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


def load_matrix(filepath, matrix_name) -> pd.DataFrame:
    """
    Read a weight matrix from a file.

    Args:
        filepath (str): Path to the file containing the weight matrix.

    Returns:
        pd.DataFrame: DataFrame containing the weight matrix.
    """
    print(f"Loading matrix: {matrix_name} from {filepath}")
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
    Compare all pairs of matrices and compute distance metrics.

    Args:
        matrices (list): List of tuples (name, filepath) for matrices.

    Returns:
        list: List of dictionaries with comparison results.
    """
    results = []
    for i in range(len(matrices)):
        for j in range(i + 1, len(matrices)):
            name1, filepath1 = matrices[i]
            name2, filepath2 = matrices[j]
            metrics = compute_distance_metrics(filepath1, filepath2, name1, name2)
            results.append({"matrix1": name1, "matrix2": name2, **metrics})
    return results


def main():
    # Get input files from snakemake
    base_matrix_files = snakemake.input.base_matrices
    predicted_matrix_files = snakemake.input.predicted_matrices
    combinations_file = snakemake.input.combinations
    output_file = snakemake.output.results

    # Load combinations to get parameter info for each predicted matrix
    combinations_df = pd.read_csv(combinations_file, sep="\t")

    # Prepare all matrices for comparison
    matrices = []

    # Add base matrices
    for matrix_file in base_matrix_files:
        name = os.path.splitext(os.path.basename(matrix_file))[0]
        matrices.append((name, matrix_file))

    # Add predicted matrices with parameter information in the name
    for matrix_file in predicted_matrix_files:
        # Extract combination ID from filename
        filename = os.path.basename(matrix_file)
        combination_id = int(filename.split("_")[-1].split(".")[0])

        # Get parameter info for this combination
        combo_row = combinations_df[
            combinations_df["combination_id"] == combination_id
        ].iloc[0]

        # Create informative name including parameters
        if combo_row["score_transform"] == "threshold":
            name = f"predicted_comb{combination_id}_steps{combo_row['steps']}_threshold{combo_row['threshold']}"
        else:  # sigmoid
            name = f"predicted_comb{combination_id}_steps{combo_row['steps']}_sigmoid_steep{combo_row['steepness']}_mid{combo_row['midpoint']}"

        matrices.append((name, matrix_file))

    print(f"Comparing {len(matrices)} matrices: {[name for name, _ in matrices]}")
    results = compare_matrices(matrices)

    # Create output DataFrame
    results_df = pd.DataFrame(results)

    # Reorder columns to match requested format
    results_df = results_df[
        [
            "matrix1",
            "matrix2",
            "euclidean_distance",
            "wasserstein_distance",
            "energy_distance",
        ]
    ]

    # Save results
    results_df.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    main()
