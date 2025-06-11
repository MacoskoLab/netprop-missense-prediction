rule download_am_scores:
    output:
        score=f"results/{run}/alphamissense/am_score.tsv"
    params:
        gene=lambda wildcards, config: config['perturb_params']['gene']
        variant=lambda wildcards, config: config['perturb_params']['variant']
    conda:
        "workflow/envs/download_am_scores.yml"
    script:
        "workflow/scripts/download_am_scores.py"
