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


rule compute_all_weight_matrices_distances:
    """
    Evaluate all weight matrices against each other using multiple distance metrics:
    - Euclidean distance
    - Wasserstein distance  
    - Energy distance
    
    Takes a list of weight matrix files and performs pairwise comparison,
    outputting results in TSV format.
    """
    input:
        base_matrices=[
            f"results/{run}/perturbation/real_unperturbed_weights.tsv",
            f"results/{run}/perturbation/real_perturbed_weights.tsv",
        ],
        predicted_matrices=expand(
            f"results/{run}/perturbation/predicted_perturbed_weights_{{combination_id}}.h5",
            combination_id=get_combination_ids(),
        ),
        combinations=f"results/{run}/perturbation/grid_combinations.tsv",
    output:
        results=f"results/{run}/evaluation/weight_matrices_comparison.tsv",
    message:
        "Evaluating weight matrices distances for all parameter combinations"
    conda:
        f"{ENVS_DIR}/evaluation.yml"
    script:
        f"{SCRIPTS_DIR}/evaluation/compute_weight_matrices_distances.py"


rule plot_weight_matrices_distances:
    """
    Create a visualization of weight matrices distances.
    Takes the TSV output from evaluate_all_weight_matrices and creates
    a bar plot showing the three distance metrics for each matrix comparison.
    Outputs both JPEG (static) and HTML (interactive) formats.
    """
    input:
        distances=f"results/{run}/evaluation/weight_matrices_comparison.tsv",
    output:
        expand(
            f"results/{run}/evaluation/weight_matrices_distances_plot.{{ext}}",
            ext=["jpeg", "html"],
        ),
    message:
        "Plotting weight matrices distances for all parameter combinations"
    conda:
        f"{ENVS_DIR}/plotting.yml"
    script:
        f"{SCRIPTS_DIR}/evaluation/plot_weight_matrices_distances.py"
