rule perturb:
    input:
        genie3_links=f"results/{run}/genie3/output.tsv",
        am_scores=f"results/{run}/get_am_scores/am_scores.tsv",
        perturbations_list=lambda wildcards: config['perturb_params']['perturbations_list']
    output:
        perturb_scores=f"results/{run}/perturb/scores.tsv",
    params:
        steps=lambda wildcards: config['perturb_params']['steps'],
        steepness=lambda wildcards: config['perturb_params']['steepness'],
        midpoint=lambda wildcards: config['perturb_params']['midpoint'],
        threshold=lambda wildcards: config['perturb_params']['threshold'],
        pathogenicity_score_transform_method=lambda wildcards: config['perturb_params']['pathogenicity_score_transform_method']
    conda:
        "../envs/perturb.yml"
    log:
        f"results/{run}/perturb/perturb.log"
    script:
        "../scripts/perturb.py"
