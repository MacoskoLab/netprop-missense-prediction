rule preprocess_expr_data:
    """
    Preprocess single-cell expression data using scanpy.
    This rule processes the full gene expression matrix, applies quality control and normalization,
    then filters cells by genotype (MT or WT, excluding NG) and saves processed expression matrices.
    Uses a Jupyter notebook for transparent and reproducible preprocessing.
    """
    input:
        expr=config["input_data"]["gene_expression_matrix"],
        metadata=config["input_data"]["cell_metadata"],
        perturbations_list=lambda wildcards: config["perturbation_algorithm"][
            "perturbations_list"
        ],
    output:
        real_unperturbed_processed_expr=f"results/{run}/preprocessing/real_unperturbed_preprocessed_expr.tsv",
        real_perturbed_processed_expr=f"results/{run}/preprocessing/real_perturbed_preprocessed_expr.tsv",
    message:
        "Preprocessing single-cell expression data"
    conda:
        f"{ENVS_DIR}/preprocess.yml"
    script:
        f"{SCRIPTS_DIR}/preprocess.py"
