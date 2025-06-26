rule download_am_scores:
    """Download AM scores from the specified source."""
    input:
        variants=lambda wildcards: config["alphamissense_api"]["variants"],
    output:
        am_scores=f"results/{run}/perturbation/am_scores.tsv",
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
    output:
        perturbed_weights=f"results/{run}/perturbation/predicted_perturbed_weights.h5",
    conda:
        f"{ENVS_DIR}/perturb.yml"
    script:
        f"{SCRIPTS_DIR}/perturbation/perturb_algo.py"
