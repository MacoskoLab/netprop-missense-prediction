rule get_am_scores:
    input:
        perturbations=lambda wildcards: config['perturb_params']['perturbations'],
    output:
        am_scores=f"results/{run}/get_am_scores/am_scores.tsv",
    conda:
        "../envs/get_am_scores.yml"
    log:
        f"results/{run}/get_am_scores/get_am_scores.log"
    script:
        "../scripts/get_am_scores.py"
