from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.script import Snakemake

    snakemake: Snakemake
    snakemake = None  # type: ignore

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import seaborn as sns
from plotly.subplots import make_subplots
from scipy.stats import energy_distance, wasserstein_distance

warnings.filterwarnings("ignore")

# Set style
plt.style.use("default")
sns.set_palette("husl")


def create_output_directory(output_dir: str):
    """Create output directory if it doesn't exist."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)


def plot_edge_weight_difference_scatter(
    real_diffs: np.ndarray, pred_diffs: np.ndarray, metrics: dict, output_path: str
):
    """Create scatter plot of real vs predicted edge weight differences."""

    fig, ax = plt.subplots(figsize=(10, 8))

    # Create scatter plot
    scatter = ax.scatter(real_diffs, pred_diffs, alpha=0.6, s=20)

    # Add diagonal line (perfect correlation)
    min_val = min(np.min(real_diffs), np.min(pred_diffs))
    max_val = max(np.max(real_diffs), np.max(pred_diffs))
    ax.plot(
        [min_val, max_val],
        [min_val, max_val],
        "r--",
        alpha=0.8,
        linewidth=2,
        label="Perfect correlation",
    )

    # Add labels and title
    ax.set_xlabel("Real Edge Weight Differences", fontsize=12)
    ax.set_ylabel("Predicted Edge Weight Differences", fontsize=12)
    ax.set_title(
        "Real vs Predicted Edge Weight Differences", fontsize=14, fontweight="bold"
    )

    # Add correlation info as text
    pearson_r = metrics.get("pearson_r", np.nan)
    spearman_r = metrics.get("spearman_r", np.nan)
    n_edges = metrics.get("n_edges", 0)

    text_str = f"n_edges = {n_edges}\n"
    text_str += f"Pearson r = {pearson_r:.3f}\n"
    text_str += f"Spearman r = {spearman_r:.3f}"

    ax.text(
        0.05,
        0.95,
        text_str,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

    # Add legend
    ax.legend()

    # Set equal aspect ratio
    ax.set_aspect("equal", adjustable="box")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_edge_weight_difference_distributions(
    real_diffs: np.ndarray, pred_diffs: np.ndarray, output_path: str
):
    """Create histogram comparing distributions of edge weight differences."""

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

    # Plot real differences
    ax1.hist(
        real_diffs, bins=50, alpha=0.7, color="blue", edgecolor="black", linewidth=0.5
    )
    ax1.set_title(
        "Distribution of Real Edge Weight Differences", fontsize=14, fontweight="bold"
    )
    ax1.set_xlabel("Edge Weight Difference", fontsize=12)
    ax1.set_ylabel("Frequency", fontsize=12)
    ax1.axvline(
        np.mean(real_diffs),
        color="red",
        linestyle="--",
        linewidth=2,
        label=f"Mean = {np.mean(real_diffs):.3f}",
    )
    ax1.axvline(
        np.median(real_diffs),
        color="green",
        linestyle="--",
        linewidth=2,
        label=f"Median = {np.median(real_diffs):.3f}",
    )
    ax1.legend()

    # Plot predicted differences
    ax2.hist(
        pred_diffs, bins=50, alpha=0.7, color="orange", edgecolor="black", linewidth=0.5
    )
    ax2.set_title(
        "Distribution of Predicted Edge Weight Differences",
        fontsize=14,
        fontweight="bold",
    )
    ax2.set_xlabel("Edge Weight Difference", fontsize=12)
    ax2.set_ylabel("Frequency", fontsize=12)
    ax2.axvline(
        np.mean(pred_diffs),
        color="red",
        linestyle="--",
        linewidth=2,
        label=f"Mean = {np.mean(pred_diffs):.3f}",
    )
    ax2.axvline(
        np.median(pred_diffs),
        color="green",
        linestyle="--",
        linewidth=2,
        label=f"Median = {np.median(pred_diffs):.3f}",
    )
    ax2.legend()

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_comparison_metrics_bar(metrics: dict, output_path: str):
    """Create bar plot of comparison metrics."""

    # Select metrics to plot
    metrics_to_plot = {
        "Pearson r": metrics.get("pearson_r", np.nan),
        "Spearman r": metrics.get("spearman_r", np.nan),
        "MSE": metrics.get("mse", np.nan),
        "MAE": metrics.get("mae", np.nan),
        "Cosine Dist": metrics.get("cosine_distance", np.nan),
    }

    # Remove NaN values
    metrics_to_plot = {k: v for k, v in metrics_to_plot.items() if not np.isnan(v)}

    if not metrics_to_plot:
        print("No valid metrics to plot")
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    names = list(metrics_to_plot.keys())
    values = list(metrics_to_plot.values())

    bars = ax.bar(
        names,
        values,
        color=["skyblue", "lightcoral", "lightgreen", "lightsalmon", "plum"],
    )

    # Add value labels on bars
    for bar, value in zip(bars, values):
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            height,
            f"{value:.3f}",
            ha="center",
            va="bottom",
            fontsize=10,
        )

    ax.set_title(
        "Edge Weight Difference Comparison Metrics", fontsize=14, fontweight="bold"
    )
    ax.set_ylabel("Metric Value", fontsize=12)

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_edge_weight_change_comparison(
    real_df: pd.DataFrame, pred_df: pd.DataFrame, output_path: str
):
    """Create side-by-side comparison of edge weight changes."""

    # Merge dataframes to get common edges
    merged = pd.merge(
        real_df[["regulatoryGene", "targetGene", "weight_diff"]],
        pred_df[["regulatoryGene", "targetGene", "weight_diff"]],
        on=["regulatoryGene", "targetGene"],
        suffixes=("_real", "_pred"),
        how="inner",
    )

    if len(merged) == 0:
        print("No common edges found for comparison plot")
        return

    # Sort by absolute real difference to highlight top changes
    merged["abs_real_diff"] = np.abs(merged["weight_diff_real"])
    top_edges = merged.nlargest(20, "abs_real_diff")

    # Create edge labels
    top_edges["edge_label"] = (
        top_edges["regulatoryGene"] + " → " + top_edges["targetGene"]
    )

    fig, ax = plt.subplots(figsize=(14, 8))

    x_pos = np.arange(len(top_edges))
    width = 0.35

    bars1 = ax.bar(
        x_pos - width / 2,
        top_edges["weight_diff_real"],
        width,
        label="Real",
        color="skyblue",
        alpha=0.8,
    )
    bars2 = ax.bar(
        x_pos + width / 2,
        top_edges["weight_diff_pred"],
        width,
        label="Predicted",
        color="lightcoral",
        alpha=0.8,
    )

    ax.set_xlabel("Edges (Regulator → Target)", fontsize=12)
    ax.set_ylabel("Edge Weight Difference", fontsize=12)
    ax.set_title(
        "Top 20 Edge Weight Changes: Real vs Predicted", fontsize=14, fontweight="bold"
    )
    ax.set_xticks(x_pos)
    ax.set_xticklabels(top_edges["edge_label"], rotation=45, ha="right", fontsize=9)
    ax.legend()

    # Add horizontal line at y=0
    ax.axhline(y=0, color="black", linestyle="-", alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_edge_weight_difference_plotly(
    real_diffs: np.ndarray, pred_diffs: np.ndarray, metrics: dict, output_path: str
):
    """Create plotly scatter plot of real vs predicted edge weight differences in the style of plot_network_comparison.py."""

    # Create single scatter plot
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=real_diffs,
            y=pred_diffs,
            mode="markers",
            opacity=0.7,
            showlegend=False,
        )
    )

    # Add identity line from -1 to 1
    fig.add_trace(
        go.Scatter(
            x=[-1, 1],
            y=[-1, 1],
            mode="lines",
            line=dict(dash="dash", color="black"),
            showlegend=False,
        )
    )

    # Update layout
    fig.update_layout(
        title_text="Real vs. Predicted Edge Weight Differences",
        title_font_size=20,
        title_x=0.5,
        height=500,
        width=600,
        margin=dict(t=60),
    )

    # Update axes with fixed range from -1 to 1
    fig.update_xaxes(
        title_text="Real edge weight difference",
        title_font_size=16,
        tickfont_size=14,
        range=[-1, 1],
        showgrid=False,
    )
    fig.update_yaxes(
        title_text="Predicted edge weight difference",
        title_font_size=16,
        tickfont_size=14,
        range=[-1, 1],
        showgrid=False,
    )

    # Save plot
    fig.write_image(output_path, scale=2)


def compute_edge_weight_differences_for_barplot(
    real_control_df: pd.DataFrame,
    real_perturbed_df: pd.DataFrame,
    pred_perturbed_df: pd.DataFrame,
) -> dict:
    """Compute edge weight differences for all three network comparisons."""

    # Compute real edge weight differences (experimental - control)
    real_control_vs_exp = pd.merge(
        real_control_df[["regulatoryGene", "targetGene", "weight"]],
        real_perturbed_df[["regulatoryGene", "targetGene", "weight"]],
        on=["regulatoryGene", "targetGene"],
        suffixes=("_control", "_exp"),
        how="inner",
    )
    real_exp_diffs = (
        real_control_vs_exp["weight_exp"] - real_control_vs_exp["weight_control"]
    ).values

    # Compute predicted vs experimental differences
    pred_vs_exp = pd.merge(
        pred_perturbed_df[["regulatoryGene", "targetGene", "weight"]],
        real_perturbed_df[["regulatoryGene", "targetGene", "weight"]],
        on=["regulatoryGene", "targetGene"],
        suffixes=("_pred", "_exp"),
        how="inner",
    )
    pred_exp_diffs = (pred_vs_exp["weight_pred"] - pred_vs_exp["weight_exp"]).values

    # Compute predicted vs control differences
    pred_vs_control = pd.merge(
        pred_perturbed_df[["regulatoryGene", "targetGene", "weight"]],
        real_control_df[["regulatoryGene", "targetGene", "weight"]],
        on=["regulatoryGene", "targetGene"],
        suffixes=("_pred", "_control"),
        how="inner",
    )
    pred_control_diffs = (
        pred_vs_control["weight_pred"] - pred_vs_control["weight_control"]
    ).values

    # Remove NaN values
    real_exp_diffs = real_exp_diffs[~np.isnan(real_exp_diffs)]
    pred_exp_diffs = pred_exp_diffs[~np.isnan(pred_exp_diffs)]
    pred_control_diffs = pred_control_diffs[~np.isnan(pred_control_diffs)]

    # Compute distance metrics for each comparison
    # For "control vs experimental", we compare the actual edge weight differences
    comparisons = {
        "control_vs_experimental": {
            "wasserstein": wasserstein_distance(
                np.zeros_like(real_exp_diffs), real_exp_diffs
            ),
            "euclidean": np.linalg.norm(real_exp_diffs),
            "e_distance": energy_distance(
                np.zeros_like(real_exp_diffs), real_exp_diffs
            ),
        },
        "predicted_vs_experimental": {
            "wasserstein": wasserstein_distance(real_exp_diffs, pred_exp_diffs),
            "euclidean": np.linalg.norm(real_exp_diffs - pred_exp_diffs),
            "e_distance": energy_distance(real_exp_diffs, pred_exp_diffs),
        },
        "predicted_vs_control": {
            "wasserstein": wasserstein_distance(
                np.zeros_like(pred_control_diffs), pred_control_diffs
            ),
            "euclidean": np.linalg.norm(pred_control_diffs),
            "e_distance": energy_distance(
                np.zeros_like(pred_control_diffs), pred_control_diffs
            ),
        },
    }

    return comparisons


def create_edge_weight_difference_barplot(
    real_control_df: pd.DataFrame,
    real_perturbed_df: pd.DataFrame,
    pred_perturbed_df: pd.DataFrame,
    output_path: str,
):
    """Create bar plots for edge weight difference distance metrics in the style of plot_network_comparison.py."""

    # Compute comparisons
    comparisons = compute_edge_weight_differences_for_barplot(
        real_control_df, real_perturbed_df, pred_perturbed_df
    )

    # Create data structure similar to plot_network_comparison.py
    results_data = [
        {
            "network_1": "control",
            "network_2": "experimental perturbation",
            "wasserstein": comparisons["control_vs_experimental"]["wasserstein"],
            "euclidean": comparisons["control_vs_experimental"]["euclidean"],
            "e_distance": comparisons["control_vs_experimental"]["e_distance"],
        },
        {
            "network_1": "predicted perturbation",
            "network_2": "experimental perturbation",
            "wasserstein": comparisons["predicted_vs_experimental"]["wasserstein"],
            "euclidean": comparisons["predicted_vs_experimental"]["euclidean"],
            "e_distance": comparisons["predicted_vs_experimental"]["e_distance"],
        },
        {
            "network_1": "predicted perturbation",
            "network_2": "control",
            "wasserstein": comparisons["predicted_vs_control"]["wasserstein"],
            "euclidean": comparisons["predicted_vs_control"]["euclidean"],
            "e_distance": comparisons["predicted_vs_control"]["e_distance"],
        },
    ]

    # Create bar plots - 3 plots stacked vertically without subplot titles
    metrics = ["wasserstein", "euclidean", "e_distance"]
    metric_titles = ["Wasserstein", "Euclidean", "Energy"]

    # Color scheme to match the original plot
    colors = ["#6366f1", "#ef4444", "#10b981"]  # Blue, Red, Green

    fig_bars = make_subplots(rows=3, cols=1, vertical_spacing=0.08, shared_xaxes=True)

    for i, (metric, title) in enumerate(zip(metrics, metric_titles), 1):
        # Get values for this metric across all comparisons
        values = [row[metric] for row in results_data]

        # Create more semantic labels
        semantic_labels = [
            "Experimental Change<br>(Control vs. Experimental)",
            "Prediction Error<br>(Predicted vs. Experimental)",
            "Predicted Change<br>(Control vs. Predicted)",
        ]

        fig_bars.add_trace(
            go.Bar(
                x=semantic_labels,
                y=values,
                name=title,
                showlegend=False,
                text=[f"{val:.3f}" for val in values],
                textposition="outside",
                textfont=dict(size=18, color="black"),
                marker=dict(color=colors[i - 1]),
            ),
            row=i,
            col=1,
        )

    # Update layout for bar plots
    fig_bars.update_layout(
        title_text="Aggregate Edge Weight Distances",
        height=900,
        width=1000,
        title_font_size=26,
        title_y=0.98,
        font=dict(size=16),
        margin=dict(t=80, b=80, l=100, r=80),
        plot_bgcolor="rgba(240,240,240,0.2)",  # Light gray background
        paper_bgcolor="white",
    )

    # Update y-axis labels to be the distance type and remove gridlines
    for i, (metric, title) in enumerate(zip(metrics, metric_titles), 1):
        # Get max value for this metric to set appropriate range
        values = [row[metric] for row in results_data]
        max_val = max(values) if max(values) > 0 else 1.0

        fig_bars.update_yaxes(
            title_text=title,
            title_font_size=20,
            tickfont_size=16,
            showgrid=False,
            zeroline=False,
            range=[0, max_val * 1.15],  # Add 15% padding at top for text labels
            row=i,
            col=1,
        )

    # Update x-axes - only show labels on bottom subplot, remove gridlines and ticks
    for i in range(1, 4):
        if i < 3:  # Top and middle subplots
            fig_bars.update_xaxes(
                showticklabels=False,
                showgrid=False,
                zeroline=False,
                ticks="",
                row=i,
                col=1,
            )
        else:  # Bottom subplot
            fig_bars.update_xaxes(
                title_font_size=20,
                tickfont_size=20,
                showgrid=False,
                zeroline=False,
                row=i,
                col=1,
            )

    # Save bar plot
    fig_bars.write_image(output_path, scale=2)


def main():
    """Main function to create all plots for edge weight differences comparison."""

    # Create output directory
    output_dir = snakemake.output.plots_dir
    create_output_directory(output_dir)

    # Load data
    print("Loading edge weight difference data...")
    real_diff_df = pd.read_csv(snakemake.input.real_edge_diffs, sep="\t")
    pred_diff_df = pd.read_csv(snakemake.input.pred_edge_diffs, sep="\t")
    comparison_results = pd.read_csv(snakemake.input.comparison_results, sep="\t")

    # Load original network files for barplot
    real_control_df = pd.read_csv(snakemake.input.real_control_weights, sep="\t")
    real_perturbed_df = pd.read_csv(snakemake.input.real_perturbed_weights, sep="\t")
    pred_perturbed_df = pd.read_csv(snakemake.input.pred_perturbed_weights, sep="\t")

    # Get metrics
    metrics = comparison_results.iloc[0].to_dict()

    # Merge data for plotting
    merged = pd.merge(
        real_diff_df[["regulatoryGene", "targetGene", "weight_diff"]],
        pred_diff_df[["regulatoryGene", "targetGene", "weight_diff"]],
        on=["regulatoryGene", "targetGene"],
        suffixes=("_real", "_pred"),
        how="inner",
    )

    if len(merged) == 0:
        print("No common edges found for plotting")
        return

    # Remove NaN values
    valid_mask = ~(
        np.isnan(merged["weight_diff_real"]) | np.isnan(merged["weight_diff_pred"])
    )
    merged_clean = merged[valid_mask]

    if len(merged_clean) == 0:
        print("No valid edge differences found for plotting")
        return

    real_diffs = merged_clean["weight_diff_real"].values
    pred_diffs = merged_clean["weight_diff_pred"].values

    print(f"Creating plots for {len(real_diffs)} edge differences...")

    # Create plots
    plot_edge_weight_difference_scatter(
        real_diffs,
        pred_diffs,
        metrics,
        f"{output_dir}/edge_weight_differences_scatter.png",
    )

    plot_edge_weight_difference_distributions(
        real_diffs,
        pred_diffs,
        f"{output_dir}/edge_weight_differences_distributions.png",
    )

    plot_comparison_metrics_bar(metrics, f"{output_dir}/comparison_metrics.png")

    plot_edge_weight_change_comparison(
        real_diff_df, pred_diff_df, f"{output_dir}/top_edge_weight_changes.png"
    )

    plot_edge_weight_difference_plotly(
        real_diffs,
        pred_diffs,
        metrics,
        f"{output_dir}/edge_weight_differences_plotly.png",
    )

    # Create the distance metrics barplot in the style of plot_network_comparison.py
    create_edge_weight_difference_barplot(
        real_control_df,
        real_perturbed_df,
        pred_perturbed_df,
        f"{output_dir}/edge_weight_difference_distance_barplot.png",
    )

    print(f"All plots saved to: {output_dir}")


if __name__ == "__main__":
    main()
