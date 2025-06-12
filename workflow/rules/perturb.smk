rule perturb:
    input:
        genie3_links=f"results/{run}/genie3/output.tsv",
        am_scores=f"results/{run}/get_am_scores/am_scores.tsv",
    output:
        perturb_scores=f"results/{run}/perturb/scores.tsv",
    params:
        k=lambda wildcards: config['perturb_params']['k'],
        t=lambda wildcards: config['perturb_params']['t'],
        pathogenicity_score_transform=lambda wildcards: config['perturb_params']['pathogenicity_score_transform']
    conda:
        "../envs/perturb.yml"
    log:
        f"results/{run}/perturb/perturb.log"
    script:
        "../scripts/perturb.py"
