"""
Utility functions for handling parameter combinations in the perturbation algorithm.
"""

import itertools


def get_combination_ids(config):
    """
    Get all combination IDs from the config.

    Args:
        config: Snakemake config dictionary

    Returns:
        list: List of combination IDs (integers)
    """
    perturb_config = config["perturbation_algorithm"]

    combination_id = 0
    combination_ids = []

    steps_list = perturb_config["steps"]
    transform_methods = perturb_config["score_transform"]

    for transform_method in transform_methods:
        for steps in steps_list:
            if transform_method == "threshold":
                threshold_values = perturb_config["threshold_params"]["threshold"]
                for threshold in threshold_values:
                    combination_ids.append(combination_id)
                    combination_id += 1
            elif transform_method == "sigmoid":
                steepness_values = perturb_config["sigmoid_params"]["steepness"]
                midpoint_values = perturb_config["sigmoid_params"]["midpoint"]
                for steepness, midpoint in itertools.product(
                    steepness_values, midpoint_values
                ):
                    combination_ids.append(combination_id)
                    combination_id += 1

    return combination_ids


def get_total_combinations_count(config):
    """
    Get the total number of combinations.

    Args:
        config: Snakemake config dictionary

    Returns:
        int: Total number of combinations
    """
    return len(get_combination_ids(config))
