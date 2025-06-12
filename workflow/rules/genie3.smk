rule genie3:
    threads: workflow.cores
    input:
        expr=config['genie3_params']['scrnaseq_input']
    output:
        links=f"results/{run}/genie3/output.tsv"
    params:
        n_trees=lambda wildcards: config['genie3_params']['n_trees'],
    conda:
        "../envs/genie3.yml"
    log:
        f"results/{run}/genie3/genie3.log"
    script:
        "../scripts/genie3.R"

