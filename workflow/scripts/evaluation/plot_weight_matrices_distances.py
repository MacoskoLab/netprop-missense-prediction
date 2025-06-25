"""
Plot weight matrices distances visualization.
Creates a bar plot showing distance metrics for matrix comparisons.
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.script import Snakemake

    snakemake: Snakemake
    snakemake = None  # type: ignore

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from utils import plotting_utils


def create_comparison_labels(df):
    """Create readable labels for matrix comparisons."""
    labels = []
    for _, row in df.iterrows():
        matrix1 = (
            row["matrix1"]
            .replace("_weights", "")
            .replace("real_", "")
            .replace("predicted_", "")
        )
        matrix2 = (
            row["matrix2"]
            .replace("_weights", "")
            .replace("real_", "")
            .replace("predicted_", "")
        )

        # Create more readable labels
        label_map = {
            "unperturbed": "Control",
            "perturbed": "Experimental",
            "predicted_perturbed": "Predicted",
        }

        label1 = label_map.get(matrix1, matrix1.title())
        label2 = label_map.get(matrix2, matrix2.title())

        labels.append(f"{label1} vs. {label2}")

    return labels


def plot_distances(df, jpeg_path, html_path):
    """Create a bar plot of distance metrics using Plotly."""
    # Create comparison labels
    comparison_labels = create_comparison_labels(df)

    # Define colors for each metric
    colors = ["#6c5ce7", "#e74c3c", "#00b894"]  # Purple, Red, Green
    metrics = ["wasserstein_distance", "euclidean_distance", "energy_distance"]
    metric_labels = ["Wasserstein", "Euclidean", "Energy"]

    # Create subplots
    fig = make_subplots(
        rows=3,
        cols=1,
        subplot_titles=metric_labels,
        vertical_spacing=0.12,
        shared_xaxes=True,
    )

    # Add traces for each metric
    for i, (metric, color, label) in enumerate(zip(metrics, colors, metric_labels)):
        values = df[metric].tolist()

        # Add bar trace
        fig.add_trace(
            go.Bar(
                x=comparison_labels,
                y=values,
                name=label,
                marker_color=color,
                text=[f"{v:.3f}" for v in values],
                textposition="outside",
                textfont=dict(size=12, color="black"),
                showlegend=False,
            ),
            row=i + 1,
            col=1,
        )

        # Update y-axis for each subplot
        fig.update_yaxes(
            title_text=label,
            title_font=dict(size=12, color="black"),
            range=[0, max(values) * 1.15],
            row=i + 1,
            col=1,
        )

        # Only show x-axis labels on the bottom subplot
        if i < 2:
            fig.update_xaxes(showticklabels=False, row=i + 1, col=1)
        else:
            fig.update_xaxes(tickangle=0, title_font=dict(size=12), row=i + 1, col=1)

    # Update layout
    fig.update_layout(
        title=dict(text="Edge Weight Difference Distance", x=0.5, font=dict(size=16)),
        height=800,
        width=800,
        showlegend=False,
    )

    # The plotting_utils template handles most styling, just customize what's needed
    fig.update_xaxes(showgrid=False, zeroline=False)
    fig.update_yaxes(showgrid=True, gridcolor="lightgray", zeroline=False)

    # Save as JPEG (requires kaleido)
    fig.write_image(jpeg_path, width=800, height=800, scale=2, format="jpeg")

    # Save as HTML (interactive)
    fig.write_html(html_path)


def main():
    # Get input and output files from snakemake
    distances_file = snakemake.input.distances
    output_files = snakemake.output

    # Extract JPEG and HTML paths from the output list
    jpeg_output = [f for f in output_files if f.endswith(".jpeg")][0]
    html_output = [f for f in output_files if f.endswith(".html")][0]

    # Load the distances data
    df = pd.read_csv(distances_file, sep="\t")

    # Create the plot
    plot_distances(df, jpeg_output, html_output)


if __name__ == "__main__":
    main()
