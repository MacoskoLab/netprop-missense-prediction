rule create_ground_truth_unperturbed_network:
    threads: workflow.cores
    input:
        expr=config['genie3_params']['scrnaseq_input']
    output:
        links=f"results/{run}/real_unperturbed_network.tsv"
    params:
        n_trees=lambda wildcards: config['genie3_params']['n_trees'],
    conda:
        "../envs/genie3.yml"
    script:
        "../scripts/genie3.R"


rule create_ground_truth_perturbed_network:
    threads: workflow.cores
    input:
        expr=config['compare_networks_params']['ground_truth_scrnaseq']
    output:
        links=f"results/{run}/real_perturbed_network.tsv"
    params:
        n_trees=lambda wildcards: config['genie3_params']['n_trees'],
    conda:
        "../envs/genie3.yml"
    script:
        "../scripts/genie3.R"
