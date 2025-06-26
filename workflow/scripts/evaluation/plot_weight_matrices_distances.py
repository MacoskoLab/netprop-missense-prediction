"""
Plot weight matrices distances visualization.
Creates a bar plot showing distance metrics for matrix comparisons.
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.script import Snakemake

    snakemake: Snakemake
    snakemake = None  # type: ignore

import os
import sys

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Add the scripts directory to Python path for utils import
scripts_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, scripts_dir)

from utils.plotting_utils import *


def create_comparison_labels(df):
    """Create readable labels for matrix comparisons."""
    labels = []
    for index, row in df.iterrows():
        print(
            f"Processing row {index}: matrix1={row['matrix1']}, matrix2={row['matrix2']}"
        )

        # Create labels based on original matrix names
        def get_label(matrix_name):
            if "predicted_perturbed" in matrix_name:
                return "Predicted"
            elif "real_unperturbed" in matrix_name:
                return "Control"
            elif "real_perturbed" in matrix_name:
                return "Experimental"
            else:
                return matrix_name.title()

        label1 = get_label(row["matrix1"])
        label2 = get_label(row["matrix2"])

        print(f"  Mapped labels: label1={label1}, label2={label2}")

        final_label = f"{label1} vs. {label2}"
        print(f"  Final label: {final_label}")
        labels.append(final_label)

    print(f"All labels: {labels}")
    return labels


def plot_distances(df, jpeg_path, html_path):
    """Create a bar plot of distance metrics using Plotly."""
    # Debug: Print dataframe info
    print(f"DataFrame shape: {df.shape}")
    print(f"DataFrame columns: {list(df.columns)}")
    print("DataFrame contents:")
    print(df)

    # Create comparison labels
    comparison_labels = create_comparison_labels(df)
    print(f"Comparison labels: {comparison_labels}")

    # Define colors for each comparison
    colors = ["#6c5ce7", "#e74c3c", "#00b894"]  # Purple, Red, Green
    metrics = ["wasserstein_distance", "euclidean_distance", "energy_distance"]
    metric_labels = ["Wasserstein", "Euclidean", "Energy"]

    # Create subplots
    fig = make_subplots(
        rows=3,
        cols=1,
        vertical_spacing=0.1,
        shared_xaxes=False,
    )

    # Add traces for each metric
    for i, (metric, metric_label) in enumerate(zip(metrics, metric_labels)):
        values = df[metric].tolist()

        # Assign colors based on position in the list
        bar_colors = [colors[j % len(colors)] for j in range(len(values))]

        # Format text labels and determine if we need scientific notation
        use_scientific = any(v < 0.001 and v > 0 for v in values if v != 0)
        text_labels = []
        for v in values:
            if v == 0:
                text_labels.append("0")
            elif use_scientific and v < 0.001:
                text_labels.append(f"{v:.2e}")
            elif v < 0.01:
                text_labels.append(f"{v:.4f}")
            else:
                text_labels.append(f"{v:.3f}")

        print(
            f"Metric {metric}: values={values}, use_scientific={use_scientific}, text_labels={text_labels}"
        )

        # Add bar trace for this metric
        fig.add_trace(
            go.Bar(
                x=comparison_labels,
                y=values,
                name=metric_label,
                marker_color=bar_colors,  # Use different colors for each bar
                text=text_labels,
                textposition="outside",
                textfont=dict(size=11, color="black"),
                showlegend=False,
            ),
            row=i + 1,
            col=1,
        )

        # Update y-axis for each subplot with appropriate formatting
        max_val = max(values) if values else 1
        if use_scientific:
            tick_format = ".2e"
        elif max_val < 0.01:
            tick_format = ".4f"
        else:
            tick_format = ".3f"

        fig.update_yaxes(
            title_text=metric_label,
            title_font=dict(size=14, color="black"),
            range=[0, max_val * 1.15] if max_val > 0 else [0, 1],
            tickformat=tick_format,
            row=i + 1,
            col=1,
        )

        # Update x-axis for each subplot
        fig.update_xaxes(
            title_font=dict(size=12),
            row=i + 1,
            col=1,
        )

    # Update layout
    fig.update_layout(
        title=dict(text="Edge Weight Difference Distance", x=0.5, font=dict(size=16)),
        height=800,
        width=800,
        showlegend=False,
        bargap=0.1,  # Gap between bars
    )

    # Update grid styling
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

    print(f"Reading distances file: {distances_file}")

    # Extract JPEG and HTML paths from the output list
    jpeg_output = [f for f in output_files if f.endswith(".jpeg")][0]
    html_output = [f for f in output_files if f.endswith(".html")][0]

    # Load the distances data
    df = pd.read_csv(distances_file, sep="\t")

    print(f"Loaded dataframe with shape: {df.shape}")
    print("Raw dataframe:")
    print(df)
    print("\nDataframe dtypes:")
    print(df.dtypes)

    # Create the plot
    plot_distances(df, jpeg_output, html_output)


if __name__ == "__main__":
    main()
