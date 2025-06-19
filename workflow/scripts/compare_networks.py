from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.script import Snakemake

    snakemake: Snakemake
    snakemake = None  # type: ignore


import numpy as np
import pandas as pd
from scipy.stats import energy_distance, wasserstein_distance


# Load networks from Snakemake inputs
def load_networks():
    """Load network data from Snakemake inputs."""
    real_unperturbed = pd.read_csv(snakemake.input.real_unperturbed, sep="\t")
    real_perturbed = pd.read_csv(snakemake.input.real_perturbed, sep="\t")
    predicted_perturbed = pd.read_csv(snakemake.input.pred, sep="\t")

    return real_unperturbed, real_perturbed, predicted_perturbed


def compute_distances(df1, df2, name1, name2):
    """Compute distances between two networks."""
    gene_pairs1 = set(
        df1[["regulatoryGene", "targetGene"]].itertuples(index=False, name=None)
    )
    gene_pairs2 = set(
        df2[["regulatoryGene", "targetGene"]].itertuples(index=False, name=None)
    )

    if gene_pairs1 != gene_pairs2:
        print(f"Some edges do not match between networks in {name1} and {name2}.")

    # Merge on regulatoryGene and targetGene
    merged = pd.merge(
        df1,
        df2,
        on=["regulatoryGene", "targetGene"],
        how="outer",
        suffixes=("_1", "_2"),
    )

    print(len(merged), "edges in the merged network.")
    na_edges = merged[merged["weight_1"].isna() | merged["weight_2"].isna()]
    if not na_edges.empty:
        print("Edges with missing weights:")
        print(na_edges[["regulatoryGene", "targetGene", "weight_1", "weight_2"]])
    else:
        print("No edges have missing weights.")

    merged[["weight_1", "weight_2"]] = merged[["weight_1", "weight_2"]].fillna(0)

    # Extract weight vectors as numpy arrays
    w1 = merged["weight_1"].astype(float).to_numpy()
    w2 = merged["weight_2"].astype(float).to_numpy()

    # Compute distances
    euclidean = np.linalg.norm(w2 - w1)
    wass = wasserstein_distance(w1, w2)
    e_distance = energy_distance(w1, w2)

    return {
        "euclidean": euclidean,
        "wasserstein": wass,
        "e_distance": e_distance,
        "weights_1": w1,
        "weights_2": w2,
    }


def perform_comparisons(real_unperturbed, real_perturbed, predicted_perturbed):
    """Perform all three pairwise network comparisons."""
    comparisons = {}
    comparisons["real_unperturbed_vs_real_perturbed"] = compute_distances(
        real_unperturbed, real_perturbed, "real unperturbed", "real perturbed"
    )
    comparisons["real_unperturbed_vs_predicted_perturbed"] = compute_distances(
        real_unperturbed, predicted_perturbed, "real unperturbed", "predicted perturbed"
    )
    comparisons["real_perturbed_vs_predicted_perturbed"] = compute_distances(
        real_perturbed, predicted_perturbed, "real perturbed", "predicted perturbed"
    )

    return comparisons


def create_results_dataframe(comparisons):
    """Create results DataFrame from comparison results."""
    comparison_names = [
        "real_unperturbed_vs_real_perturbed",
        "real_unperturbed_vs_predicted_perturbed",
        "real_perturbed_vs_predicted_perturbed",
    ]

    network_pairs = [
        ("real unperturbed", "real perturbed"),
        ("real unperturbed", "predicted perturbed"),
        ("real perturbed", "predicted perturbed"),
    ]

    results_data = []
    for name, (net1, net2) in zip(comparison_names, network_pairs):
        comp = comparisons[name]
        results_data.append(
            {
                "network_1": net1,
                "network_2": net2,
                "wasserstein": comp["wasserstein"],
                "euclidean": comp["euclidean"],
                "e_distance": comp["e_distance"],
            }
        )

    return pd.DataFrame(results_data)


def save_results(results_df, output_path):
    """Save results DataFrame to TSV file."""
    results_df.to_csv(output_path, sep="\t", index=False)
    print(f"Results saved to: {output_path}")


def main():
    """Main function to orchestrate the network comparison workflow."""
    # Load network data
    real_unperturbed, real_perturbed, predicted_perturbed = load_networks()

    # Perform all pairwise comparisons
    comparisons = perform_comparisons(
        real_unperturbed, real_perturbed, predicted_perturbed
    )

    # Create results DataFrame
    results_df = create_results_dataframe(comparisons)

    # Save results
    save_results(results_df, snakemake.output.comp)


if __name__ == "__main__":
    main()
