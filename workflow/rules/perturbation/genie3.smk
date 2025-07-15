rule build_genie3_weights:
    """
    Build the weight matrix for either unperturbed or perturbed expression data
    depending on the `{state}` wildcard (allowed values: ``unperturbed``, ``perturbed``).
    """
    threads: 1
    input:
        expr=f"results/{run}/preprocessing/real_{{state}}_preprocessed_expr.tsv",
    output:
        weights=f"results/{run}/perturbation/real_{{state}}_weights.tsv",
    message:
        f"Building GENIE3 weight matrix for {{wildcards.state}} expression data."
    benchmark:
        f"results/{run}/perturbation/benchmark_genie3_real_{{state}}_weights.txt"
    conda:
        f"{ENVS_DIR}/genie3.yml"
    script:
        f"{SCRIPTS_DIR}/perturbation/genie3_weights.R"
