
rule run_genie3:
    """
    Build either the unperturbed or perturbed ground truth network
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
        links=f"results/{run}/real_{{state}}_network.tsv",
    params:
        tree_method=lambda wc: config["genie3_params"]["tree_method"],
        K=lambda wc: config["genie3_params"]["K"],
        n_trees=lambda wc: config["genie3_params"]["n_trees"],
    conda:
        "../envs/genie3.yml"
    script:
        "../scripts/genie3.R"
