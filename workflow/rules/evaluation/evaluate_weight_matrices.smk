rule evaluate_all_weight_matrices:
    """
    Evaluate all weight matrices against each other using multiple distance metrics:
    - Euclidean distance
    - Wasserstein distance  
    - Energy distance
    
    Takes a list of weight matrix files and performs pairwise comparison,
    outputting results in TSV format.
    """
    input:
        matrices=[
            f"results/{run}/perturbation/real_unperturbed_weights.tsv",
            f"results/{run}/perturbation/real_perturbed_weights.tsv",
            f"results/{run}/perturbation/predicted_perturbed_weights.tsv",
        ],
    output:
        results=f"results/{run}/evaluation/weight_matrices_comparison.tsv",
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
    conda:
        f"{ENVS_DIR}/evaluation.yml"
    script:
        f"{SCRIPTS_DIR}/evaluation/plot_weight_matrices_distances.py"
