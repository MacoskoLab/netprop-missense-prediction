"""
Generate grid combinations for perturbation algorithm parameters.
This script creates a TSV file mapping combination IDs to parameter sets,
ensuring no redundant combinations between threshold and sigmoid methods.
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.script import Snakemake

    snakemake: Snakemake
    snakemake = None  # type: ignore

import itertools

import pandas as pd


def generate_combinations(config):
    """Generate all valid parameter combinations."""
    perturb_config = config["perturbation_algorithm"]

    combinations = []
    combination_id = 0

    # Get common parameters
    steps_list = perturb_config["steps"]
    transform_methods = perturb_config["score_transform"]

    for transform_method in transform_methods:
        for steps in steps_list:
            if transform_method == "threshold":
                # Generate combinations for threshold method
                threshold_values = perturb_config["threshold_params"]["threshold"]
                for threshold in threshold_values:
                    combination = {
                        "combination_id": combination_id,
                        "steps": steps,
                        "score_transform": transform_method,
                        "threshold": threshold,
                        "steepness": None,  # Not used for threshold
                        "midpoint": None,  # Not used for threshold
                    }
                    combinations.append(combination)
                    combination_id += 1

            elif transform_method == "sigmoid":
                # Generate combinations for sigmoid method
                steepness_values = perturb_config["sigmoid_params"]["steepness"]
                midpoint_values = perturb_config["sigmoid_params"]["midpoint"]

                for steepness, midpoint in itertools.product(
                    steepness_values, midpoint_values
                ):
                    combination = {
                        "combination_id": combination_id,
                        "steps": steps,
                        "score_transform": transform_method,
                        "threshold": None,  # Not used for sigmoid
                        "steepness": steepness,
                        "midpoint": midpoint,
                    }
                    combinations.append(combination)
                    combination_id += 1

    return combinations


def main():
    """Main function to generate grid combinations."""
    config = snakemake.config
    combinations = generate_combinations(config)

    # Convert to DataFrame
    df = pd.DataFrame(combinations)

    # Save to TSV
    output_path = snakemake.output.combinations
    df.to_csv(output_path, sep="\t", index=False)

    print(f"Generated {len(combinations)} parameter combinations")
    print(f"Grid combinations saved to: {output_path}")

    # Print summary
    threshold_count = len(
        [c for c in combinations if c["score_transform"] == "threshold"]
    )
    sigmoid_count = len([c for c in combinations if c["score_transform"] == "sigmoid"])

    print(f"Threshold combinations: {threshold_count}")
    print(f"Sigmoid combinations: {sigmoid_count}")


if __name__ == "__main__":
    main()
