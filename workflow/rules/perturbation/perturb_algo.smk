import sys
import os

# Add the scripts directory to Python path for utils import
scripts_dir = os.path.join(workflow.basedir, "scripts")
sys.path.insert(0, scripts_dir)

from utils.combination_utils import get_combination_ids


rule generate_grid_combinations:
    """Generate parameter combinations for perturbation algorithm grid search."""
    output:
        combinations=f"results/{run}/perturbation/grid_combinations.tsv",
    message:
        "Generating parameter combinations for perturbation algorithm grid search"
    script:
        f"{SCRIPTS_DIR}/perturbation/generate_grid_combinations.py"


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
        f"{SCRIPTS_DIR}/perturbation/download_am_scores.py"


rule simple_perturb_algo:
    """Run the simple perturbation algorithm using AM scores."""
    input:
        genie3_weights=f"results/{run}/perturbation/real_unperturbed_weights.tsv",
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
        f"{SCRIPTS_DIR}/perturbation/perturb_algo.py"
