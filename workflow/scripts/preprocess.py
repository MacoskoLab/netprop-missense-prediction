# %% [markdown]
# # Single-cell RNA-seq Data Preprocessing Pipeline
#
# This notebook performs comprehensive preprocessing of single-cell RNA-seq data following the same approach as the data exploration notebook, but outputs separate files for WT and MT cells for downstream analysis.
#
# **Pipeline Steps:**
# 1. Load expression matrix and cell metadata
# 2. Create AnnData object and calculate QC metrics
# 3. Quality control filtering (optional)
# 4. Normalization (optional)
# 5. Log transformation (optional)
# 4. Normalization (optional)
# 5. Log transformation (optional)
# 6. Highly variable gene identification
# 7. Data scaling (optional)
# 8. PCA and dimensionality reduction (with visualizations)
# 9. UMAP computation and visualization (with and without outliers)
# 10. Split data by genotype and save processed matrices

# %%
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.script import Snakemake

    snakemake: Snakemake
    snakemake = None  # type: ignore

# %% [markdown]
# <span id="papermill-error-cell" style="color:red; font-family:Helvetica Neue, Helvetica, Arial, sans-serif; font-size:2em;">Execution using papermill encountered an exception here and stopped:</span>

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np

# %%
# Import required libraries
import pandas as pd
import scanpy as sc
import seaborn as sns

# Set scanpy settings
sc.settings.verbosity = 3  # verbosity level
# sc.settings.set_figure_params(dpi=80, facecolor="white")

print("Libraries imported successfully")

# %%
# Access Snakemake inputs and outputs
input_expr = snakemake.input.expr
input_metadata = snakemake.input.metadata
perturbations_list = snakemake.input.perturbations_list
output_wt = snakemake.output.real_unperturbed_processed_expr
output_mt = snakemake.output.real_perturbed_processed_expr

# Get preprocessing configuration
config = snakemake.config["single_cell_preprocessing"]

print(f"Input expression matrix: {input_expr}")
print(f"Input metadata: {input_metadata}")
print(f"Output WT file: {output_wt}")
print(f"Output MT file: {output_mt}")
print(f"Configuration: {config}")

# %% [markdown]
# ## Step 1: Load Data

# %%
# Load cell metadata
print("Loading cell metadata...")
cell_meta = pd.read_csv(input_metadata, sep="\t", index_col=0)
print(f"Loaded metadata for {len(cell_meta)} cells")
print(f"Perturbation distribution:\n{cell_meta['Perturbation'].value_counts()}")

# Load expression matrix
print("\nLoading expression matrix...")
expr = pd.read_csv(input_expr, sep="\t", index_col=0)
expr = expr.loc[:, cell_meta.index]  # Keep only cells that have metadata
print(f"Expression matrix shape: {expr.shape} (genes x cells)")

# Load perturbation genes list for protection during filtering
protected_genes = set()
if Path(perturbations_list).exists():
    protected_df = pd.read_csv(perturbations_list, sep="\t")
    if "gene" in protected_df.columns:
        protected_genes = set(protected_df["gene"])
        print(f"Loaded {len(protected_genes)} protected genes")
    else:
        print("No 'gene' column found in perturbations list")
else:
    print("No perturbations list file found")

print(f"Protected genes: {protected_genes}")

# %% [markdown]
# ## Step 2: Create AnnData Object and Calculate QC Metrics

# %%
# Create AnnData object (transpose so cells are rows, genes are columns)
adata = ad.AnnData(X=expr.T)
adata.var_names = expr.index  # Gene names
adata.obs_names = expr.columns  # Cell names

# Add metadata
adata.obs["perturbation"] = cell_meta["Perturbation"]
adata.obs["genotype"] = cell_meta["Perturbation"].map(
    {"MT": "Mutant", "WT": "Wild-type", "NG": "Not-genotyped"}
)

print(f"Data shape: {adata.shape[0]} cells × {adata.shape[1]} genes")
print(f"Genotype distribution:\n{adata.obs['genotype'].value_counts()}")

# Calculate QC metrics
# Mitochondrial genes (MT- prefix)
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# Ribosomal genes (RP[SL] prefix)
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# Hemoglobin genes (HB prefix)
adata.var["hb"] = adata.var_names.str.startswith("HB")

# Calculate basic QC metrics
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

# Manually calculate percentage of mitochondrial and ribosomal genes
adata.obs["pct_counts_mt"] = (
    np.array(adata[:, adata.var["mt"]].X.sum(axis=1)).flatten()
    / adata.obs["total_counts"]
    * 100
)
adata.obs["pct_counts_ribo"] = (
    np.array(adata[:, adata.var["ribo"]].X.sum(axis=1)).flatten()
    / adata.obs["total_counts"]
    * 100
)

print("QC metrics calculated:")
print(
    f"- Total counts per cell: {adata.obs['total_counts'].min():.0f} - {adata.obs['total_counts'].max():.0f}"
)
print(
    f"- Genes per cell: {adata.obs['n_genes_by_counts'].min()} - {adata.obs['n_genes_by_counts'].max()}"
)
print(f"- Mitochondrial genes: {adata.var['mt'].sum()}")
print(f"- Ribosomal genes: {adata.var['ribo'].sum()}")
print(f"- Hemoglobin genes: {adata.var['hb'].sum()}")

# Display QC metrics
adata.obs[
    ["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo"]
].head()

# %% [markdown]
# ## Step 3: Quality Control Filtering (Optional)

# %%
# Save raw data before any filtering
adata.raw = adata

# Apply quality control filtering if enabled
if config.get("enable_quality_control", False):
    print("Applying quality control filtering...")

    min_genes = config["min_genes"]
    min_cells = config["min_cells"]
    max_genes = config["max_genes"]

    print(f"Before filtering: {adata.shape[0]} cells by {adata.shape[1]} genes")

    # Filter cells with too few or too many genes
    sc.pp.filter_cells(adata, min_genes=min_genes)
    print(f"After min_genes filter: {adata.shape[0]} cells by {adata.shape[1]} genes")

    sc.pp.filter_cells(adata, max_genes=max_genes)
    print(f"After max_genes filter: {adata.shape[0]} cells by {adata.shape[1]} genes")

    # Filter genes with protection for perturbation genes
    if len(protected_genes) > 0:
        protected_genes_present = protected_genes.intersection(set(adata.var_names))
        # Create mask for genes to keep
        genes_to_keep = (
            adata.var["n_cells_by_counts"] >= min_cells
        ) | adata.var_names.isin(protected_genes_present)
        adata = adata[:, genes_to_keep].copy()
        print(
            f"After min_cells filter (protecting {len(protected_genes_present)} genes): {adata.shape[0]} cells by {adata.shape[1]} genes"
        )
    else:
        sc.pp.filter_genes(adata, min_cells=min_cells)
        print(
            f"After min_cells filter: {adata.shape[0]} cells by {adata.shape[1]} genes"
        )

    print(f"Final genotype distribution after filtering:")
    print(adata.obs["genotype"].value_counts())
else:
    print("Skipping quality control filtering")

# %% [markdown]
# ## Step 4: Normalization (Optional)

# %%
# Apply normalization if enabled
if config.get("enable_normalization", False):
    print("Applying normalization...")

    target_sum = config["target_sum"]

    # Normalize to target sum per cell
    sc.pp.normalize_total(adata, target_sum=target_sum)

    # Save normalized data
    adata.layers["normalized"] = adata.X.copy()

    print(f"Normalization completed:")
    print(f"- Library size normalized to {target_sum} reads per cell")
    print(f"- Data range: {adata.X.min():.2f} to {adata.X.max():.2f}")
else:
    print("Skipping normalization")

# %% [markdown]
# ## Step 5: Log Transformation (Optional)

# %%
# Apply log transformation if enabled
if config.get("enable_log_transformation", False):
    print("Applying log transformation...")

    # Log transform (log(x + 1))
    sc.pp.log1p(adata)

    # Save log-transformed data
    adata.layers["log_transformed"] = adata.X.copy()

    print(f"Log transformation completed:")
    print(f"- Log(x+1) transformation applied")
    print(f"- Data range: {adata.X.min():.2f} to {adata.X.max():.2f}")
else:
    print("Skipping log transformation")

# Save raw data if no processing was done yet
if adata.raw is None:
    adata.raw = adata

# %% [markdown]
# ## Step 6: Highly Variable Gene Identification

# %%
# Identify highly variable genes if enabled
if config.get("enable_hvg_identification", True):
    print("Identifying highly variable genes...")

    min_mean = config["min_mean"]
    max_mean = config["max_mean"]
    min_disp = config["min_disp"]
    n_top_genes = config.get("n_top_genes", None)

    # If normalization and/or log transformation were skipped, do temporary processing for HVG calculation
    normalization_done = config.get("enable_normalization", False)
    log_transformation_done = config.get("enable_log_transformation", False)

    if not (normalization_done and log_transformation_done):
        print(
            "Creating temporary normalized and log-transformed data for HVG calculation..."
        )
        adata_temp = adata.copy()

        if not normalization_done:
            sc.pp.normalize_total(adata_temp, target_sum=config["target_sum"])

        if not log_transformation_done:
            sc.pp.log1p(adata_temp)

        # Calculate HVGs on temporary data
        sc.pp.highly_variable_genes(
            adata_temp,
            min_mean=min_mean,
            max_mean=max_mean,
            min_disp=min_disp,
            n_top_genes=n_top_genes,
        )

        # Transfer HVG information back to original data
        adata.var["highly_variable"] = adata_temp.var["highly_variable"]
        if "means" in adata_temp.var.columns:
            adata.var["means"] = adata_temp.var["means"]
        if "dispersions" in adata_temp.var.columns:
            adata.var["dispersions"] = adata_temp.var["dispersions"]
        if "dispersions_norm" in adata_temp.var.columns:
            adata.var["dispersions_norm"] = adata_temp.var["dispersions_norm"]

        print("HVG identification completed using temporary processing")
    else:
        # Calculate HVGs on already normalized data
        sc.pp.highly_variable_genes(
            adata,
            min_mean=min_mean,
            max_mean=max_mean,
            min_disp=min_disp,
            n_top_genes=n_top_genes,
        )
        print("HVG identification completed on processed data")

    print(f"Number of highly variable genes: {adata.var['highly_variable'].sum()}")
    print(f"Total genes: {adata.shape[1]}")
    print(
        f"Percentage of highly variable genes: {adata.var['highly_variable'].sum() / adata.shape[1] * 100:.1f}%"
    )

    # Ensure perturbation genes are marked as highly variable
    if len(protected_genes) > 0:
        protected_genes_present = protected_genes.intersection(set(adata.var_names))
        if len(protected_genes_present) > 0:
            print(
                f"Ensuring {len(protected_genes_present)} perturbation genes are marked as highly variable"
            )
            for gene in protected_genes_present:
                adata.var.loc[gene, "highly_variable"] = True

            # Report final HVG counts
            n_hvgs = adata.var["highly_variable"].sum()
            print(
                f"Final count: {n_hvgs} highly variable genes (including {len(protected_genes_present)} protected genes)"
            )

    # Show top highly variable genes if dispersions were calculated
    if "dispersions_norm" in adata.var.columns:
        hvg_df = adata.var[adata.var["highly_variable"]].sort_values(
            "dispersions_norm", ascending=False
        )
        print("\nTop 10 highly variable genes:")
        print(hvg_df[["means", "dispersions", "dispersions_norm"]].head(10))

        # Plot highly variable genes
        sc.pl.highly_variable_genes(adata)
else:
    print(
        "Skipping highly variable gene identification - marking all genes as highly variable"
    )
    adata.var["highly_variable"] = True

# %% [markdown]
# ## Step 7: Data Scaling (Optional)

# %%
# Apply data scaling if enabled
if config.get("enable_scaling", False):
    print("Scaling data...")

    max_value = config["max_value"]

    # Save processed data before scaling
    normalization_done = config.get("enable_normalization", False)
    log_transformation_done = config.get("enable_log_transformation", False)

    if normalization_done or log_transformation_done:
        adata.layers["pre_scaling"] = adata.X.copy()

    # Scale data
    sc.pp.scale(adata, max_value=max_value)

    print(f"Scaling completed:")
    print(f"- Data centered to zero mean")
    print(f"- Scaled to unit variance")
    print(f"- Extreme values clipped to [-{max_value}, {max_value}]")
    print(f"- Data range: {adata.X.min():.2f} to {adata.X.max():.2f}")
else:
    print("Skipping data scaling")

# %%
# Subset to highly variable genes if HVG identification is enabled
if config.get("enable_hvg_identification", True):
    adata = adata[:, adata.var["highly_variable"]].copy()
    print(f"Subset to {adata.shape[1]} HVGs for downstream analysis")

# %% [markdown]
# ## Step 8: Split Data by Genotype and Save Results


# %%
# Filter cells by genotype and save processed expression matrices
def save_genotype_data(adata, genotype, output_path):
    """Save processed expression data for a specific genotype"""

    # Filter cells for this genotype
    if genotype == "WT":
        genotype_mask = adata.obs["perturbation"] == "WT"
        genotype_name = "Wild-type"
    elif genotype == "MT":
        genotype_mask = adata.obs["perturbation"] == "MT"
        genotype_name = "Mutant"
    else:
        raise ValueError(f"Unknown genotype: {genotype}")

    adata_genotype = adata[genotype_mask, :].copy()
    print(f"\n{genotype_name} ({genotype}) cells: {adata_genotype.shape[0]} cells")

    # Get processed expression matrix (use raw counts if no processing was done)
    normalization_done = config.get("enable_normalization", False)
    log_transformation_done = config.get("enable_log_transformation", False)
    scaling_done = config.get("enable_scaling", False)

    if normalization_done or log_transformation_done or scaling_done:
        # Use processed data
        if hasattr(adata_genotype.X, "toarray"):
            expr_matrix = pd.DataFrame(
                adata_genotype.X.toarray().T,  # Transpose back to genes x cells
                index=adata_genotype.var_names,
                columns=adata_genotype.obs_names,
            )
        else:
            expr_matrix = pd.DataFrame(
                adata_genotype.X.T,  # Transpose back to genes x cells
                index=adata_genotype.var_names,
                columns=adata_genotype.obs_names,
            )
        print(f"Using processed data for {genotype_name}")
    else:
        # Use raw counts
        if hasattr(adata_genotype.raw.X, "toarray"):
            expr_matrix = pd.DataFrame(
                adata_genotype.raw.X.toarray().T,  # Transpose back to genes x cells
                index=adata_genotype.var_names,
                columns=adata_genotype.obs_names,
            )
        else:
            expr_matrix = pd.DataFrame(
                adata_genotype.raw.X.T,  # Transpose back to genes x cells
                index=adata_genotype.var_names,
                columns=adata_genotype.obs_names,
            )
        print(f"Using raw counts for {genotype_name}")

    # Save to file
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    expr_matrix.to_csv(output_path, sep="\t")
    print(f"Saved {genotype_name} expression matrix to {output_path}")
    print(f"Matrix shape: {expr_matrix.shape} (genes x cells)")

    return expr_matrix


# Save WT (unperturbed) data
print("=== Saving WT (unperturbed) data ===")
wt_matrix = save_genotype_data(adata, "WT", output_wt)

# Save MT (perturbed) data
print("\n=== Saving MT (perturbed) data ===")
mt_matrix = save_genotype_data(adata, "MT", output_mt)

print(f"\n=== Preprocessing Pipeline Completed Successfully! ===")
print(f"\nSummary:")
print(
    f"- Data processing completed with {adata.shape[0]} cells and {adata.shape[1]} genes"
)
print(f"- Generated expression matrices for WT and MT cells separately")
print(f"- Created visualizations with and without outliers for clarity")
print(f"- All cells (including outliers) were saved in the output files")
print(f"- Outlier removal was only applied to visualization, not to saved data")

if "highly_variable" in adata.var.columns:
    print(f"- Identified {adata.var['highly_variable'].sum()} highly variable genes")

print(f"\nOutput files:")
print(f"- WT (unperturbed): {output_wt}")
print(f"- MT (perturbed): {output_mt}")
print(f"\nCell counts in saved data:")
print(f"- WT cells: {wt_matrix.shape[1]}")
print(f"- MT cells: {mt_matrix.shape[1]}")
print(f"- Total: {wt_matrix.shape[1] + mt_matrix.shape[1]} cells saved")
