# Validate genie3_params configuration
def validate_genie3_params():
    """Validate that threshold and top_n_links are mutually exclusive"""
    genie3_params = config.get("genie3_params", {})
    threshold = genie3_params.get("threshold")
    top_n_links = genie3_params.get("top_n_links")

    # Check if both are non-null
    if (
        threshold is not None
        and threshold != "null"
        and top_n_links is not None
        and top_n_links != "null"
    ):
        raise ValueError(
            "Error in config: genie3_params 'threshold' and 'top_n_links' are mutually exclusive. "
            "Please set only one of them to a value and the other to null."
        )

    print("GENIE3 parameters validation passed", flush=True)


validate_genie3_params()


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
        weights=f"results/{run}/real_{{state}}_weights.tsv",
    params:
        tree_method=lambda wc: config["genie3_params"]["tree_method"],
        K=lambda wc: config["genie3_params"]["K"],
        n_trees=lambda wc: config["genie3_params"]["n_trees"],
    conda:
        "../envs/genie3.yml"
    script:
        "../scripts/genie3_weights.R"


rule get_genie3_links:
    """
    Generate the ranked list of regulatory links from the weight matrix
    for either unperturbed or perturbed networks.
    """
    input:
        weights=f"results/{run}/real_{{state}}_weights.tsv",
    output:
        links=f"results/{run}/real_{{state}}_network.tsv",
    params:
        threshold=lambda wc: config["genie3_params"]["threshold"],
        top_n_links=lambda wc: config["genie3_params"]["top_n_links"],
    conda:
        "../envs/genie3.yml"
    script:
        "../scripts/genie3_links.R"
