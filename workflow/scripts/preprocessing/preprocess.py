"""
Comprehensive single-cell preprocessing pipeline using scanpy.

This script performs:
- Quality control filtering
- Normalization and log-transformation
- Highly variable gene identification
- Scaling and PCA
- Cell clustering
- Pseudobulk aggregation
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.script import Snakemake

    snakemake: Snakemake
    snakemake = None  # type: ignore

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


def load_expression_matrix(filepath):
    """
    Load expression matrix from file.

    Parameters:
    -----------
    filepath : str
        Path to the expression matrix file

    Returns:
    --------
    pd.DataFrame
        Expression matrix with genes as rows and cells as columns
    """
    print(f"Loading expression matrix from {filepath}")
    return pd.read_csv(filepath, sep="\t", index_col=0, header=0)


def quality_control_filtering(
    adata, min_genes=200, min_cells=3, max_genes=5000, mt_threshold=20
):
    """
    Filter out low-quality cells and genes with low expression.

    Parameters:
    -----------
    adata : anndata.AnnData
        Input AnnData object
    min_genes : int
        Minimum number of genes expressed per cell
    min_cells : int
        Minimum number of cells expressing each gene
    max_genes : int
        Maximum number of genes per cell (filter potential doublets)
    mt_threshold : float
        Maximum percentage of mitochondrial genes per cell (unused)

    Returns:
    --------
    anndata.AnnData
        Quality-filtered AnnData object
    """
    print("Performing quality control filtering...")

    # Calculate basic QC metrics
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

    print(f"Initial data shape: {adata.shape}")

    # Filter cells with too few or too many genes
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, max_genes=max_genes)

    # Filter genes expressed in too few cells
    sc.pp.filter_genes(adata, min_cells=min_cells)

    print(f"After QC filtering: {adata.shape}")
    return adata


def normalize_and_log_transform(adata, target_sum=1e4):
    """
    Normalize counts per cell and log-transform.

    Parameters:
    -----------
    adata : anndata.AnnData
        Input AnnData object
    target_sum : float
        Target sum for normalization (counts per cell)

    Returns:
    --------
    anndata.AnnData
        Normalized and log-transformed AnnData object
    """
    print("Normalizing and log-transforming data...")

    # Save raw counts
    adata.raw = adata

    # Normalize to target sum per cell
    sc.pp.normalize_total(adata, target_sum=target_sum)

    # Log-transform
    sc.pp.log1p(adata)

    print("Normalization and log-transformation completed")
    return adata


def identify_highly_variable_genes(
    adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=None
):
    """
    Identify highly variable genes on normalized, log-transformed data.

    Parameters:
    -----------
    adata : anndata.AnnData
        Normalized, log-transformed AnnData object
    min_mean : float
        Minimum mean expression for HVG selection
    max_mean : float
        Maximum mean expression for HVG selection
    min_disp : float
        Minimum normalized dispersion for HVG selection
    n_top_genes : int or None
        Number of top variable genes to select

    Returns:
    --------
    anndata.AnnData
        AnnData object with HVG information
    """
    print("Identifying highly variable genes...")

    # Find highly variable genes
    sc.pp.highly_variable_genes(
        adata,
        min_mean=min_mean,
        max_mean=max_mean,
        min_disp=min_disp,
        n_top_genes=n_top_genes,
    )

    n_hvgs = adata.var["highly_variable"].sum()
    print(
        f"Identified {n_hvgs} highly variable genes out of {adata.n_vars} total genes"
    )

    return adata


def scale_data(adata, max_value=10):
    """
    Scale data for highly variable genes.

    Parameters:
    -----------
    adata : anndata.AnnData
        Input AnnData object with HVG information
    max_value : float
        Maximum value after scaling

    Returns:
    --------
    anndata.AnnData
        Scaled AnnData object
    """
    print("Scaling data for highly variable genes...")

    # Keep only highly variable genes for scaling
    adata_hvg = adata[:, adata.var.highly_variable].copy()

    # Scale data
    sc.pp.scale(adata_hvg, max_value=max_value)

    print(f"Scaled data shape: {adata_hvg.shape}")
    return adata_hvg


def run_pca(adata, n_comps=50):
    """
    Run PCA on scaled, variable-gene data.

    Parameters:
    -----------
    adata : anndata.AnnData
        Scaled AnnData object
    n_comps : int
        Number of principal components

    Returns:
    --------
    anndata.AnnData
        AnnData object with PCA embeddings
    """
    print(f"Running PCA with {n_comps} components...")

    # Principal component analysis
    sc.tl.pca(adata, svd_solver="arpack", n_comps=n_comps)

    print("PCA completed")
    return adata


def cluster_cells(adata, n_neighbors=10, n_pcs=40, resolution=0.5):
    """
    Cluster cells using PCA embeddings.

    Parameters:
    -----------
    adata : anndata.AnnData
        AnnData object with PCA embeddings
    n_neighbors : int
        Number of neighbors for neighborhood graph
    n_pcs : int
        Number of PCs to use for neighborhood graph
    resolution : float
        Resolution for Leiden clustering

    Returns:
    --------
    anndata.AnnData
        AnnData object with cluster assignments
    """
    print("Clustering cells...")

    # Compute neighborhood graph
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    # Leiden clustering
    sc.tl.leiden(adata, resolution=resolution)

    n_clusters = len(adata.obs["leiden"].unique())
    print(f"Identified {n_clusters} clusters")

    return adata


def create_pseudobulk(adata_raw, adata_processed, sample_col=None):
    """
    Aggregate raw counts per cluster per sample (pseudobulk).

    Parameters:
    -----------
    adata_raw : anndata.AnnData
        Raw count data
    adata_processed : anndata.AnnData
        Processed data with cluster assignments
    sample_col : str or None
        Column name for sample information

    Returns:
    --------
    pd.DataFrame
        Pseudobulk expression matrix
    """
    print("Creating pseudobulk profiles...")

    # Get cluster assignments
    clusters = adata_processed.obs["leiden"].values

    # Create sample-cluster combinations
    if sample_col and sample_col in adata_processed.obs.columns:
        samples = adata_processed.obs[sample_col].values
        pseudobulk_groups = [
            f"sample_{s}_cluster_{c}" for s, c in zip(samples, clusters)
        ]
    else:
        # If no sample info, just use clusters
        pseudobulk_groups = [f"cluster_{c}" for c in clusters]

    # Create pseudobulk matrix by summing raw counts
    pseudobulk_data = {}
    unique_groups = list(set(pseudobulk_groups))

    for group in unique_groups:
        group_mask = [g == group for g in pseudobulk_groups]
        group_counts = adata_raw.X[group_mask, :].sum(axis=0)

        # Handle sparse matrices
        if hasattr(group_counts, "A1"):
            group_counts = group_counts.A1

        pseudobulk_data[group] = group_counts

    # Create DataFrame
    pseudobulk_df = pd.DataFrame(pseudobulk_data, index=adata_raw.var_names)

    print(
        f"Created pseudobulk matrix with shape: {pseudobulk_df.shape} (genes x samples)"
    )
    return pseudobulk_df


def create_anndata_object(expr_df):
    """
    Create AnnData object from expression DataFrame.

    Parameters:
    -----------
    expr_df : pd.DataFrame
        Expression matrix with genes as rows and cells as columns

    Returns:
    --------
    anndata.AnnData
        AnnData object with proper orientation (cells x genes)
    """
    print("Creating AnnData object")

    # Transpose so that cells are observations (rows) and genes are variables (columns)
    adata = ad.AnnData(X=expr_df.T)
    adata.obs_names = expr_df.columns
    adata.var_names = expr_df.index

    # Make variable names unique
    adata.var_names_make_unique()

    print(f"Created AnnData object with shape: {adata.shape} (cells x genes)")
    return adata


def save_filtered_matrix(adata_filtered, output_path, format="tsv"):
    """
    Save filtered expression matrix to file.

    Parameters:
    -----------
    adata_filtered : anndata.AnnData
        Filtered AnnData object
    output_path : str
        Output file path
    format : str
        Output format ('tsv', 'csv', 'h5ad')
    """
    print(f"Saving filtered matrix to {output_path}")

    if format == "h5ad":
        adata_filtered.write_h5ad(output_path)
    else:
        # Convert back to DataFrame with genes as rows, cells as columns
        df_filtered = pd.DataFrame(
            adata_filtered.X.T,
            index=adata_filtered.var.index,
            columns=adata_filtered.obs.index,
        )

        df_filtered.to_hdf(output_path, key="expression", mode="w", format="fixed")

    print("Filtered matrix saved successfully")


def load_cell_metadata(filepath):
    """
    Load cell metadata file with perturbation information.

    Parameters:
    -----------
    filepath : str
        Path to the cell metadata file

    Returns:
    --------
    pd.DataFrame
        DataFrame with cell IDs and perturbation information
    """
    print(f"Loading cell metadata from {filepath}")
    metadata = pd.read_csv(filepath, sep="\t", index_col=0)
    print(f"Loaded metadata for {len(metadata)} cells")
    print(f"Perturbation types: {metadata['Perturbation'].value_counts()}")
    return metadata


def filter_cells_by_genotype(adata, metadata, genotype):
    """
    Filter AnnData object to include only cells of a specific genotype.

    Parameters:
    -----------
    adata : anndata.AnnData
        Input AnnData object
    metadata : pd.DataFrame
        Cell metadata with perturbation information
    genotype : str
        Genotype to filter for ('MT' or 'WT')

    Returns:
    --------
    anndata.AnnData
        Filtered AnnData object
    """
    print(f"Filtering cells for genotype: {genotype}")

    # Get cells that match the genotype and are present in the expression data
    genotype_cells = metadata[metadata["Perturbation"] == genotype].index
    common_cells = [cell for cell in genotype_cells if cell in adata.obs_names]

    print(f"Found {len(common_cells)} {genotype} cells in expression data")

    # Filter the AnnData object
    adata_filtered = adata[common_cells, :].copy()

    return adata_filtered


def main():
    # Get parameters from Snakemake config
    input_expr = snakemake.input.expr
    input_metadata = snakemake.input.metadata
    output_wt = snakemake.output.real_unperturbed_processed_expr
    output_mt = snakemake.output.real_perturbed_processed_expr

    # Access parameters directly from config
    preprocessing_config = snakemake.config["single_cell_preprocessing"]

    min_mean = preprocessing_config["min_mean"]
    max_mean = preprocessing_config["max_mean"]
    min_disp = preprocessing_config["min_disp"]
    n_top_genes = preprocessing_config.get("n_top_genes", None)

    # Quality control parameters
    min_genes = preprocessing_config["min_genes"]
    min_cells = preprocessing_config["min_cells"]
    max_genes = preprocessing_config["max_genes"]
    mt_threshold = preprocessing_config["mt_threshold"]

    # Normalization parameters
    target_sum = preprocessing_config["target_sum"]

    # Scaling parameters
    max_value = preprocessing_config["max_value"]

    # PCA parameters
    n_comps = preprocessing_config["n_comps"]

    # Clustering parameters
    n_neighbors = preprocessing_config["n_neighbors"]
    n_pcs = preprocessing_config["n_pcs"]
    resolution = preprocessing_config["resolution"]

    # Load expression matrix
    expr_df = load_expression_matrix(input_expr)

    # Load cell metadata
    metadata = load_cell_metadata(input_metadata)

    # Create AnnData object
    adata = create_anndata_object(expr_df)

    # Step 1: Quality control filtering (applied to all cells)
    adata = quality_control_filtering(
        adata,
        min_genes=min_genes,
        min_cells=min_cells,
        max_genes=max_genes,
        mt_threshold=mt_threshold,
    )

    # Step 2: Normalize and log-transform (applied to all cells)
    adata = normalize_and_log_transform(adata, target_sum=target_sum)

    # Step 3: Identify highly variable genes (applied to all cells)
    adata = identify_highly_variable_genes(
        adata,
        min_mean=min_mean,
        max_mean=max_mean,
        min_disp=min_disp,
        n_top_genes=n_top_genes,
    )

    # Process WT cells (unperturbed)
    print("\n=== Processing WT (unperturbed) cells ===")
    adata_wt = filter_cells_by_genotype(adata, metadata, "WT")

    # Scale, PCA, and cluster WT cells
    adata_wt_scaled = scale_data(adata_wt, max_value=max_value)
    adata_wt_scaled = run_pca(adata_wt_scaled, n_comps=n_comps)
    adata_wt_scaled = cluster_cells(
        adata_wt_scaled, n_neighbors=n_neighbors, n_pcs=n_pcs, resolution=resolution
    )

    # Create pseudobulk profiles for WT
    pseudobulk_wt = create_pseudobulk(adata_wt.raw.to_adata(), adata_wt_scaled)

    # Save WT pseudobulk matrix
    print(f"Saving WT pseudobulk matrix to {output_wt}")
    pseudobulk_wt.to_csv(output_wt, sep="\t", index=True)

    # Process MT cells (perturbed)
    print("\n=== Processing MT (perturbed) cells ===")
    adata_mt = filter_cells_by_genotype(adata, metadata, "MT")

    # Scale, PCA, and cluster MT cells
    adata_mt_scaled = scale_data(adata_mt, max_value=max_value)
    adata_mt_scaled = run_pca(adata_mt_scaled, n_comps=n_comps)
    adata_mt_scaled = cluster_cells(
        adata_mt_scaled, n_neighbors=n_neighbors, n_pcs=n_pcs, resolution=resolution
    )

    # Create pseudobulk profiles for MT
    pseudobulk_mt = create_pseudobulk(adata_mt.raw.to_adata(), adata_mt_scaled)

    # Save MT pseudobulk matrix
    print(f"Saving MT pseudobulk matrix to {output_mt}")
    pseudobulk_mt.to_csv(output_mt, sep="\t", index=True)

    print("Single-cell preprocessing pipeline completed successfully!")
    print(f"Generated WT (unperturbed) file: {output_wt}")
    print(f"Generated MT (perturbed) file: {output_mt}")


if __name__ == "__main__":
    main()
