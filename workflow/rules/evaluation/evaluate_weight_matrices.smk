def get_all_perturbation_files(wildcards):
    """Get all predicted perturbation files for evaluation."""
    import pandas as pd
    import os

    combinations_file = f"results/{wildcards.run}/perturbation/grid_combinations.tsv"

    # Check if combinations file exists (for grid search mode)
    if os.path.exists(combinations_file):
        combinations_df = pd.read_csv(combinations_file, sep="\t")
        combination_ids = combinations_df["combination_id"].tolist()
        return [
            f"results/{wildcards.run}/perturbation/predicted_perturbed_weights_{cid}.h5"
            for cid in combination_ids
        ]
    else:
        # Fallback to single file (for backward compatibility)
        return [f"results/{wildcards.run}/perturbation/predicted_perturbed_weights.h5"]


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
        predicted_matrices=get_all_perturbation_files,
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
        report(
            expand(
                f"results/{run}/evaluation/weight_matrices_distances_plot.{{ext}}",
                ext=["jpeg", "html"],
            ),
            caption="report/weight_matrices_distances.rst",
            category="Weight Matrix Evaluation",
            subcategory="Distance Metrics",
        ),
    message:
        "Plotting weight matrices distances for all parameter combinations"
    conda:
        f"{ENVS_DIR}/plotting.yml"
    script:
        f"{SCRIPTS_DIR}/evaluation/plot_weight_matrices_distances.py"
