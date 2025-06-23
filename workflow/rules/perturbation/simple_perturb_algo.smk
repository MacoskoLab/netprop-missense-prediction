rule download_am_scores:
    """Download AM scores from the specified source."""
    input:
        variants=lambda wildcards: config["am_params"]["variants"],
    output:
        am_scores=f"results/{run}/perturbation/am_scores.tsv",
    conda:
        f"../../{ENVS_DIR}/download_am_scores.yml"
    script:
        f"../../{SCRIPTS_DIR}/perturbation/download_am_scores.py"


rule run_simple_perturb_algo:
    """Run the simple perturbation algorithm using AM scores."""
    input:
        genie3_links=f"results/{run}/preprocessing/real_unperturbed_network.tsv",
        am_scores=f"results/{run}/perturbation/am_scores.tsv",
        perturbations_list=lambda wildcards: config["perturb_params"][
            "perturbations_list"
        ],
    output:
        scores=f"results/{run}/perturbation/predicted_perturbed_network.tsv",
    params:
        steps=lambda wildcards: config["perturb_params"]["steps"],
        steepness=lambda wildcards: config["perturb_params"]["steepness"],
        midpoint=lambda wildcards: config["perturb_params"]["midpoint"],
        threshold=lambda wildcards: config["perturb_params"]["threshold"],
        pathogenicity_score_transform_method=lambda wildcards: config[
            "perturb_params"
        ]["pathogenicity_score_transform_method"],
        debug=lambda wildcards: config["perturb_params"]["debug"],
    conda:
        f"../../{ENVS_DIR}/perturb.yml"
    script:
        f"../../{SCRIPTS_DIR}/perturbation/perturb_algo.py"
