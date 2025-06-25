import os
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.script import Snakemake

    snakemake: Snakemake
    snakemake = None  # type: ignore

import sys

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

sys.path.append("../utils")
from plotly_util import *


def get_network_weights(net1_name, net2_name, network_map):
    """Get weights for two networks to create scatter plots."""
    df1 = network_map[net1_name]
    df2 = network_map[net2_name]

    # Merge on regulatoryGene and targetGene
    merged = pd.merge(
        df1,
        df2,
        on=["regulatoryGene", "targetGene"],
        how="outer",
        suffixes=("_1", "_2"),
    )

    # Extract weight vectors as numpy arrays
    w1 = merged["weight_1"].astype(float).to_numpy()
    w2 = merged["weight_2"].astype(float).to_numpy()

    return w1, w2


def create_scatter_plots(results_data, network_map):
    """Create scatter plots comparing network weights."""
    # Create scatter plots - 3 subplots side by side without titles
    fig_scatter = make_subplots(
        rows=1,
        cols=3,
        horizontal_spacing=0.1,
    )

    for i, row in enumerate(results_data, 1):
        w1, w2 = get_network_weights(row["network_1"], row["network_2"], network_map)

        fig_scatter.add_trace(
            go.Scatter(
                x=w1,
                y=w2,
                mode="markers",
                opacity=0.7,
                showlegend=False,
            ),
            row=1,
            col=i,
        )

        # Add identity line starting from 0
        max_val = max(w1.max(), w2.max())
        fig_scatter.add_trace(
            go.Scatter(
                x=[0, max_val],
                y=[0, max_val],
                mode="lines",
                line=dict(dash="dash", color="black"),
                showlegend=False,
            ),
            row=1,
            col=i,
        )

    # Simple layout with big title
    fig_scatter.update_layout(
        title_text="Network Weights",
        title_font_size=20,
        title_x=0.5,
        height=450,
        width=1200,
        margin=dict(t=60),
    )

    # Update axes with labels only
    for i, row in enumerate(results_data, 1):
        x_label = row["network_1"].title()
        y_label = row["network_2"].title()
        fig_scatter.update_xaxes(
            title_text=x_label,
            title_font_size=16,
            tickfont_size=14,
            range=[0, None],
            rangemode="tozero",
            row=1,
            col=i,
        )
        fig_scatter.update_yaxes(
            title_text=y_label,
            title_font_size=16,
            tickfont_size=14,
            range=[0, None],
            rangemode="tozero",
            row=1,
            col=i,
        )

    # Save scatter plot
    scatter_path = os.path.join(
        snakemake.output.figs, "network_comparison_scatterplot.png"
    )
    fig_scatter.write_image(scatter_path, scale=2)


def create_bar_plots(results_data):
    """Create bar plots for distance metrics."""
    # Create bar plots - 3 plots stacked vertically without subplot titles
    metrics = ["wasserstein", "euclidean", "e_distance"]
    metric_titles = ["Wasserstein", "Euclidean", "Energy"]

    # Color scheme to match the edge weight difference plot
    colors = ["#6366f1", "#ef4444", "#10b981"]  # Blue, Red, Green

    fig_bars = make_subplots(rows=3, cols=1, vertical_spacing=0.08, shared_xaxes=True)

    for i, (metric, title) in enumerate(zip(metrics, metric_titles), 1):
        # Get values for this metric across all comparisons
        values = [row[metric] for row in results_data]
        labels = [
            f"{row['network_1'].replace('_', ' ').replace(' perturbation', '').title()}\nvs.\n{row['network_2'].replace('_', ' ').replace(' perturbation', '').title()}"
            for row in results_data
        ]

        fig_bars.add_trace(
            go.Bar(
                x=labels,
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
        title_text="Edge Weight Difference Distance",
        height=900,
        width=900,
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
    bar_path = os.path.join(snakemake.output.figs, "network_comparison_barplot.png")
    save_plotly_fig(fig_bars, bar_path)


def create_difference_scatter_plot(results_data, network_map):
    """Create a single scatter plot comparing network differences."""
    # Get the network data
    control = network_map["control"]
    experimental = network_map["experimental perturbation"]
    insilico = network_map["in-silico perturbation"]

    # Merge all three networks on regulatoryGene and targetGene
    merged = (
        control.merge(
            experimental,
            on=["regulatoryGene", "targetGene"],
            how="outer",
            suffixes=("_control", "_exp"),
        )
        .merge(insilico, on=["regulatoryGene", "targetGene"], how="outer")
        .rename(columns={"weight": "weight_insilico"})
    )

    # Calculate differences
    control_vs_exp_diff = merged["weight_exp"] - merged["weight_control"]
    control_vs_insilico_diff = merged["weight_insilico"] - merged["weight_control"]

    # Create single scatter plot
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=control_vs_exp_diff,
            y=control_vs_insilico_diff,
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
        title_text="Experimental vs. In-silico Perturbation",
        title_font_size=20,
        title_x=0.5,
        height=500,
        width=600,
        margin=dict(t=60),
    )

    # Update axes with fixed range from -1 to 1
    fig.update_xaxes(
        title_text="Experimental perturbation minus control",
        title_font_size=16,
        tickfont_size=14,
        range=[-1, 1],
        showgrid=False,
    )
    fig.update_yaxes(
        title_text="In-silico perturbation minus control",
        title_font_size=16,
        tickfont_size=14,
        range=[-1, 1],
        showgrid=False,
    )

    # Save difference scatter plot
    diff_path = os.path.join(
        snakemake.output.figs, "network_difference_scatterplot.png"
    )
    save_plotly_fig(fig, diff_path)


def load_data():
    """Load comparison results and network data."""
    # Load comparison results
    results_df = pd.read_csv(snakemake.input.comp, sep="\t")
    results_df[["network_1", "network_2"]] = results_df[["network_1", "network_2"]].map(
        lambda network_kind: {
            "real unperturbed": "control",
            "predicted perturbed": "in-silico perturbation",
            "real perturbed": "experimental perturbation",
        }[network_kind]
    )
    results_data = results_df.to_dict("records")

    # Load networks from Snakemake inputs to get weights for scatter plots
    real_unperturbed = pd.read_csv(snakemake.input.real_unperturbed, sep="\t")
    real_perturbed = pd.read_csv(snakemake.input.real_perturbed, sep="\t")
    predicted_perturbed = pd.read_csv(snakemake.input.pred, sep="\t")

    network_map = {
        "control": real_unperturbed,
        "experimental perturbation": real_perturbed,
        "in-silico perturbation": predicted_perturbed,
    }

    return results_data, network_map


def main():
    """Main function to orchestrate the plotting workflow."""
    # Create output directory
    os.makedirs(snakemake.output.figs, exist_ok=True)

    # Load data
    results_data, network_map = load_data()

    # Create plots
    # create_scatter_plots(results_data, network_map)
    create_bar_plots(results_data)
    create_difference_scatter_plot(results_data, network_map)


if __name__ == "__main__":
    main()
