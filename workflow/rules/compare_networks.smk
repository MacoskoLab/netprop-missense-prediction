rule compare_perturbed_networks:
    input:
        real_unperturbed=f"results/{run}/real_unperturbed_network.tsv",
        real_perturbed=f"results/{run}/real_perturbed_network.tsv",
        pred=f"results/{run}/predicted_perturbed_network.tsv",
    output:
        comp=f"results/{run}/network_comparison.tsv",
    conda:
        "../envs/compare_networks.yml"
    script:
        "../scripts/compare_networks.py"


rule plot_network_comparison:
    input:
        comp=f"results/{run}/network_comparison.tsv",
        real_unperturbed=f"results/{run}/real_unperturbed_network.tsv",
        real_perturbed=f"results/{run}/real_perturbed_network.tsv",
        pred=f"results/{run}/predicted_perturbed_network.tsv",
    output:
        figs=directory(f"results/{run}/figs/network_comparison/"),
    conda:
        "../envs/plot_networks.yml"
    script:
        "../scripts/plot_network_comparison.py"
