rule get_am_scores:
    input:
        variants=lambda wildcards: config['am_params']['variants'],
    output:
        am_scores=f"results/{run}/get_am_scores/am_scores.tsv",
    conda:
        "../envs/get_am_scores.yml"
    log:
        f"results/{run}/get_am_scores/get_am_scores.log"
    script:
        "../scripts/get_am_scores.py"
