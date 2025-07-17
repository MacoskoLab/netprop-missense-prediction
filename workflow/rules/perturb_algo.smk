rule generate_grid_combinations:
    """Generate parameter combinations for perturbation algorithm grid search."""
    output:
        combinations=f"results/{run}/perturbation/grid_combinations.tsv",
    message:
        "Generating parameter combinations for perturbation algorithm grid search"
    script:
        f"{SCRIPTS_DIR}/generate_grid_combinations.py"


def get_combination_ids():
    """Get all combination IDs from the config."""
    import itertools

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


rule download_am_scores:
    """Download AM scores from the specified source."""
    input:
        variants=lambda wildcards: config["alphamissense_api"]["variants"],
    output:
        am_scores=f"results/{run}/perturbation/am_scores.tsv",
    message:
        "Downloading AM scores"
    conda:
        f"{ENVS_DIR}/download_am_scores.yml"
    script:
        f"{SCRIPTS_DIR}/download_am_scores.py"


rule simple_perturb_algo:
    """Run the simple perturbation algorithm using AM scores."""
    input:
        genie3_weights=f"results/{run}/perturbation/real_unperturbed_weights.tsv",
        real_perturbed_weights=f"results/{run}/perturbation/real_perturbed_weights.tsv",
        am_scores=f"results/{run}/perturbation/am_scores.tsv",
        perturbations_list=lambda wildcards: config["perturbation_algorithm"][
            "perturbations_list"
        ],
        combinations=f"results/{run}/perturbation/grid_combinations.tsv",
    output:
        perturbed_weights=f"results/{run}/perturbation/predicted_perturbed_weights_{{combination_id}}.h5",
    message:
        "Running simple perturbation algorithm with combination {wildcards.combination_id}"
    params:
        combination_id="{combination_id}",
    conda:
        f"{ENVS_DIR}/perturb.yml"
    script:
        f"{SCRIPTS_DIR}/perturb_algo.py"
