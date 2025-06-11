rule perturb:
    input:
        genie3_links=f"results/{run}/genie3/output.csv",
        am_score=f"results/{run}/alphamissense/am_score.txt"
    output:
        scores=f"results/{run}/perturb/scores.tsv"
    params:
        gene=lambda wildcards, config: config['perturb_params']['gene']
        k=lambda wildcards, config: config['perturb_params']['k']
        t=lambda wildcards, config: config['perturb_params']['t']
    conda:
        "workflow/envs/perturb.yml"
    script:
        "workflow/scripts/perturb.py"
