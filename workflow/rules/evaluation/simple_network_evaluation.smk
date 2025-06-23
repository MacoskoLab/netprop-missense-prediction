rule compare_networks_simple:
    input:
        real_unperturbed=f"results/{run}/preprocessing/real_unperturbed_network.tsv",
        real_perturbed=f"results/{run}/preprocessing/real_perturbed_network.tsv",
        pred=f"results/{run}/perturbation/predicted_perturbed_network.tsv",
    output:
        comp=f"results/{run}/evaluation/network_comparison/network_comparison.tsv",
    conda:
        f"../../{ENVS_DIR}/compare_networks.yml"
    script:
        f"../../{SCRIPTS_DIR}/evaluation/compare_networks.py"


rule plot_network_comparison_simple:
    input:
        comp=f"results/{run}/evaluation/network_comparison/network_comparison.tsv",
        real_unperturbed=f"results/{run}/preprocessing/real_unperturbed_network.tsv",
        real_perturbed=f"results/{run}/preprocessing/real_perturbed_network.tsv",
        pred=f"results/{run}/perturbation/predicted_perturbed_network.tsv",
    output:
        figs=directory(f"results/{run}/evaluation/network_comparison/figs/"),
    conda:
        f"../../{ENVS_DIR}/plot_networks.yml"
    script:
        f"../../{SCRIPTS_DIR}/visualization/plot_network_comparison.py"
