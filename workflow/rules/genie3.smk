rule genie3:
    input:
        expr=config['genie3_params']['scrnaseq_input']
    output:
        links=f"results/{run}/genie3/output.tsv"
    params:
        n_trees=lambda wildcards, config: config['genie3_params']['n_trees']
        min_size=lambda wildcards, config: config['genie3_params']['min_size']
        max_depth=lambda wildcards, config: config['genie3_params']['max_depth']
    conda:
        "workflow/envs/genie3.yml"
    script:
        "workflow/scripts/genie3.r"

