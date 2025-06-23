rule compute_mock_perturbed_edge_weights_from_perturbed_expression:
    """
    Compute edge weights from expression data using the control network topology.
    This creates real perturbed and control edge weights that share the same topology.
    """
    input:
        control_network=f"results/{run}/preprocessing/real_unperturbed_network.tsv",
        perturbed_expr=lambda wc: config["genie3_params"]["real_perturbed_gene_expr"],
        control_expr=lambda wc: config["genie3_params"]["real_unperturbed_gene_expr"],
    output:
        perturbed_weights=f"results/{run}/evaluation/mock_perturbation/real_perturbed_weights_consistent_topology.tsv",
        control_weights=f"results/{run}/evaluation/mock_perturbation/real_control_weights_consistent_topology.tsv",
    conda:
        f"../../{ENVS_DIR}/compare_networks.yml"
    script:
        f"../../{SCRIPTS_DIR}/evaluation/compute_edge_weights_from_expression.py"


rule compare_edge_weight_differences_mock_perturbation:
    """
    Compare edge weight differences (Δ = w_perturbed - w_control) between
    real and predicted networks using the same topology.
    """
    input:
        real_control_weights=f"results/{run}/evaluation/mock_perturbation/real_control_weights_consistent_topology.tsv",
        real_perturbed_weights=f"results/{run}/evaluation/mock_perturbation/real_perturbed_weights_consistent_topology.tsv",
        pred_control_weights=f"results/{run}/preprocessing/real_unperturbed_network.tsv",  # Original control network
        pred_perturbed_weights=f"results/{run}/perturbation/predicted_perturbed_network.tsv",
    output:
        real_edge_diffs=f"results/{run}/evaluation/mock_perturbation/real_edge_weight_differences.tsv",
        pred_edge_diffs=f"results/{run}/evaluation/mock_perturbation/predicted_edge_weight_differences.tsv",
        comparison_results=f"results/{run}/evaluation/mock_perturbation/edge_weight_difference_comparison.tsv",
    conda:
        f"../../{ENVS_DIR}/compare_networks.yml"
    script:
        f"../../{SCRIPTS_DIR}/evaluation/compare_edge_weight_differences.py"


rule plot_edge_weight_differences_mock_perturbation:
    """
    Create plots comparing real vs predicted edge weight differences.
    """
    input:
        real_edge_diffs=f"results/{run}/evaluation/mock_perturbation/real_edge_weight_differences.tsv",
        pred_edge_diffs=f"results/{run}/evaluation/mock_perturbation/predicted_edge_weight_differences.tsv",
        comparison_results=f"results/{run}/evaluation/mock_perturbation/edge_weight_difference_comparison.tsv",
        real_control_weights=f"results/{run}/evaluation/mock_perturbation/real_control_weights_consistent_topology.tsv",
        real_perturbed_weights=f"results/{run}/evaluation/mock_perturbation/real_perturbed_weights_consistent_topology.tsv",
        pred_perturbed_weights=f"results/{run}/perturbation/predicted_perturbed_network.tsv",
    output:
        plots_dir=directory(f"results/{run}/evaluation/mock_perturbation/figs/"),
    conda:
        f"../../{ENVS_DIR}/plot_networks.yml"
    script:
        f"../../{SCRIPTS_DIR}/visualization/plot_edge_weight_differences.py"
