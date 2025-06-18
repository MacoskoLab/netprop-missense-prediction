rule download_am_scores:
    input:
        variants=lambda wildcards: config["am_params"]["variants"],
    output:
        am_scores=f"results/{run}/am_scores.tsv",
    conda:
        "../envs/download_am_scores.yml"
    script:
        "../scripts/download_am_scores.py"


rule run_perturb_algo:
    input:
        genie3_links=f"results/{run}/real_unperturbed_network.tsv",
        am_scores=f"results/{run}/am_scores.tsv",
        perturbations_list=lambda wildcards: config["perturb_params"][
            "perturbations_list"
        ],
    output:
        scores=f"results/{run}/predicted_perturbed_network.tsv",
    params:
        steps=lambda wildcards: config["perturb_params"]["steps"],
        steepness=lambda wildcards: config["perturb_params"]["steepness"],
        midpoint=lambda wildcards: config["perturb_params"]["midpoint"],
        threshold=lambda wildcards: config["perturb_params"]["threshold"],
        pathogenicity_score_transform_method=lambda wildcards: config[
            "perturb_params"
        ]["pathogenicity_score_transform_method"],
    conda:
        "../envs/perturb.yml"
    script:
        "../scripts/perturb_algo.py"
