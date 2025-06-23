rule build_genie3_weights:
    """
    Build the weight matrix for either unperturbed or perturbed expression data
    depending on the `{state}` wildcard (allowed values: ``unperturbed``, ``perturbed``).
    """
    threads: workflow.cores
    # Select expression matrix based on the requested state
    input:
        expr=lambda wc: (
            config["genie3_params"]["real_unperturbed_gene_expr"]
            if wc.state == "unperturbed"
            else config["genie3_params"]["real_perturbed_gene_expr"]
        ),
    # Wildcard `state` propagates into the output name
    output:
        weights=f"results/{run}/preprocessing/real_{{state}}_weights.tsv",
    params:
        tree_method=lambda wc: config["genie3_params"]["tree_method"],
        K=lambda wc: config["genie3_params"]["K"],
        n_trees=lambda wc: config["genie3_params"]["n_trees"],
    conda:
        f"../../{ENVS_DIR}/genie3.yml"
    script:
        f"../../{SCRIPTS_DIR}/preprocessing/genie3_weights.R"


rule get_genie3_links:
    """
    Generate the ranked list of regulatory links from the weight matrix
    for either unperturbed or perturbed networks.
    """
    input:
        weights=f"results/{run}/preprocessing/real_{{state}}_weights.tsv",
    output:
        links=f"results/{run}/preprocessing/real_{{state}}_network.tsv",
    params:
        link_filtering_method=lambda wc: config["genie3_params"][
            "link_filtering_method"
        ],
        threshold=lambda wc: config["genie3_params"]["threshold"],
        top_n_links=lambda wc: config["genie3_params"]["top_n_links"],
    conda:
        f"../../{ENVS_DIR}/genie3.yml"
    script:
        f"../../{SCRIPTS_DIR}/preprocessing/genie3_links.R"
