import sys
import os

# Add the scripts directory to Python path for utils import
scripts_dir = os.path.join(workflow.basedir, "scripts")
sys.path.insert(0, scripts_dir)

from utils.combination_utils import get_combination_ids


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
        real_perturbed_matrix=f"results/{run}/perturbation/real_perturbed_weights.tsv",
        real_unperturbed_matrix=f"results/{run}/perturbation/real_unperturbed_weights.tsv",
        matrix_combinations=f"results/{run}/perturbation/grid_combinations.tsv",
        predicted_perturbed_matrices=expand(
            f"results/{run}/perturbation/predicted_perturbed_weights_{{combination_id}}.h5",
            combination_id=get_combination_ids(config),
        ),
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
    Create visualizations of weight matrices distances.
    Creates separate plots for each parameter combination, showing the three distance 
    metrics for each matrix comparison. Outputs both JPEG (static) and HTML (interactive) 
    formats for each parameter combination.
    """
    input:
        distances=f"results/{run}/evaluation/weight_matrices_comparison.tsv",
    output:
        # Create a directory to hold all the plot files
        directory(f"results/{run}/evaluation/weight_matrices_distances_plots/"),
    message:
        "Plotting weight matrices distances for all parameter combinations"
    conda:
        f"{ENVS_DIR}/plotting.yml"
    script:
        f"{SCRIPTS_DIR}/evaluation/plot_weight_matrices_distances.py"
