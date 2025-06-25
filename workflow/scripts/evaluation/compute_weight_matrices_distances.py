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


def compute_distance_metrics(matrix1, matrix2):
    """
    Compute distance metrics between two weight matrices.

    Matrices must have exactly the same genes for both rows and columns,
    and contain no NaN values.

    Returns:
        dict: Dictionary with euclidean, wasserstein, and energy distances
    """
    # Check that matrices have exactly the same row indices (genes)
    if not matrix1.index.equals(matrix2.index):
        raise ValueError(
            f"Matrices have different row genes. "
            f"Matrix1 has {len(matrix1.index)} genes, Matrix2 has {len(matrix2.index)} genes. "
            f"Row genes must be identical."
        )

    # Check that matrices have exactly the same column indices (genes)
    if not matrix1.columns.equals(matrix2.columns):
        raise ValueError(
            f"Matrices have different column genes. "
            f"Matrix1 has {len(matrix1.columns)} genes, Matrix2 has {len(matrix2.columns)} genes. "
            f"Column genes must be identical."
        )

    # Check for NaN values in matrix1
    if matrix1.isna().any().any():
        raise ValueError("Matrix1 contains NaN values. All values must be numeric.")

    # Check for NaN values in matrix2
    if matrix2.isna().any().any():
        raise ValueError("Matrix2 contains NaN values. All values must be numeric.")

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


def load_weight_matrices(matrix_files):
    """
    Load weight matrices from a list of file paths.

    Args:
        matrix_files (list): List of file paths to weight matrices.

    Returns:
        dict: Dictionary with matrix names as keys and DataFrames as values.
    """
    matrices = {}
    for filepath in matrix_files:
        matrix_name = os.path.splitext(os.path.basename(filepath))[0]
        print(f"Loading matrix: {matrix_name} from {filepath}")
        # Only TSV inputs
        matrices[matrix_name] = pd.read_csv(filepath, sep="\t", index_col=0, header=0)
    return matrices


def compare_matrices(matrices):
    """
    Compare all pairs of matrices and compute distance metrics.

    Args:
        matrices (dict): Dictionary with matrix names as keys and DataFrames as values.

    Returns:
        list: List of dictionaries with comparison results.
    """
    results = []
    for (name1, matrix1), (name2, matrix2) in combinations(matrices.items(), 2):
        metrics = compute_distance_metrics(matrix1, matrix2)
        results.append({"matrix1": name1, "matrix2": name2, **metrics})
    return results


def main():
    # Get input files from snakemake
    matrix_files = snakemake.input.matrices
    output_file = snakemake.output.results

    # Load all matrices
    matrices = load_weight_matrices(matrix_files)
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
