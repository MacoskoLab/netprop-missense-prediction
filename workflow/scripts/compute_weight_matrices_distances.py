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
    # Assume matrices are preloaded DataFrames
    matrix1 = matrix1.fillna(0)
    matrix2 = matrix2.fillna(0)

    print(f"Comparing matrices: {name1} vs {name2}", flush=True)

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
        print(f"WARNING: {name1} contains NaN values.", flush=True)

    # Check for NaN values in matrix2
    if matrix2.isna().any().any():
        print(f"WARNING: {name2} contains NaN values.", flush=True)

    # Flatten matrices for distance computation
    m1_flat = matrix1.values.flatten()
    m2_flat = matrix2.values.flatten()

    # Compute distance metrics
    euclidean_dist = euclidean(m1_flat, m2_flat)
    print("Computed Euclidean distance", flush=True)
    wasserstein_dist = wasserstein_distance(m1_flat, m2_flat)
    print("Computed Wasserstein distance", flush=True)
    energy_dist = energy_distance(m1_flat, m2_flat)
    print("Computed Energy distance", flush=True)

    print("Computed distances", flush=True)

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
    print(f"Loading matrix: {matrix_name} from {filepath}", flush=True)
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
    # Initialize
    print("Starting weight matrix distance computation...", flush=True)
    real1_fp = snakemake.input.real_unperturbed_weights
    real2_fp = snakemake.input.real_perturbed_weights
    predicted_fps = snakemake.input.predicted_matrices
    combinations_fp = snakemake.input.combinations
    output_fp = snakemake.output.results

    # Load combinations to annotate predicted matrices
    combinations_df = pd.read_csv(combinations_fp, sep="\t")
    combo_cols = combinations_df.columns.tolist()

    # Preload all matrices into memory
    matrices = {}  # name -> DataFrame
    params = {}  # name -> dict of parameters

    # Load real matrices
    for fp in [real1_fp, real2_fp]:
        nm = os.path.splitext(os.path.basename(fp))[0]
        df = load_matrix(fp, nm).fillna(0)
        matrices[nm] = df
        params[nm] = {col: None for col in combo_cols}

    # Load predicted matrices with parameters
    for fp in predicted_fps:
        fname = os.path.basename(fp)
        comb_id = int(fname.split("_")[-1].split(".")[0])
        # find combination row
        combo_row = combinations_df.loc[
            combinations_df["combination_id"] == comb_id
        ].iloc[0]
        # construct name
        if combo_row["score_transform"] == "threshold":
            nm = f"predicted_comb{comb_id}_steps{combo_row['steps']}_threshold{combo_row['threshold']}"
        else:
            nm = (
                f"predicted_comb{comb_id}_steps{combo_row['steps']}_"
                f"sigmoid_steep{combo_row['steepness']}_mid{combo_row['midpoint']}"
            )
        df = load_matrix(fp, nm).fillna(0)
        matrices[nm] = df
        params[nm] = combo_row.to_dict()

    # Prepare comparisons: real-real and real-predicted
    real_names = [
        os.path.splitext(os.path.basename(real1_fp))[0],
        os.path.splitext(os.path.basename(real2_fp))[0],
    ]
    compare_pairs = []
    # real vs real
    compare_pairs.append((real_names[0], real_names[1]))
    # real vs predicted
    for real_nm in real_names:
        for nm in matrices:
            if nm not in real_names:
                compare_pairs.append((real_nm, nm))

    # Compute metrics for each pair
    results = []
    for m1, m2 in compare_pairs:
        metrics = compute_distance_metrics(matrices[m1], matrices[m2], m1, m2)
        row = {
            "matrix1": m1,
            "matrix2": m2,
            "wasserstein_distance": metrics["wasserstein_distance"],
            "euclidean_distance": metrics["euclidean_distance"],
            "energy_distance": metrics["energy_distance"],
        }
        # add parameter columns (from predicted matrix if present, otherwise NA)
        predicted_matrix = None
        if m1.startswith("predicted_"):
            predicted_matrix = m1
        elif m2.startswith("predicted_"):
            predicted_matrix = m2

        for col in combo_cols:
            if predicted_matrix:
                row[col] = params[predicted_matrix].get(col)
            else:
                row[col] = None
        results.append(row)

    # Create output DataFrame
    results_df = pd.DataFrame(results)

    # Reorder columns: matrix1, matrix2, wasserstein, euclidean, energy, then parameters
    all_cols = [
        "matrix1",
        "matrix2",
        "wasserstein_distance",
        "euclidean_distance",
        "energy_distance",
    ] + combo_cols
    results_df = results_df[all_cols]

    results_df[["matrix1", "matrix2"]] = results_df[["matrix1", "matrix2"]].map(
        lambda matrix_name: (
            "predicted_perturbed_weights"
            if matrix_name.startswith("predicted")
            else matrix_name
        )
    )

    # Save results
    results_df.to_csv(output_fp, sep="\t", index=False)


if __name__ == "__main__":
    main()
