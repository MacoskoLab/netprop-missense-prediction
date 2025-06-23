"""Configure Plotly templates and default settings."""

import datetime
import pathlib

import plotly.graph_objects as go
import plotly.io as pio


def save_plotly_fig(fig, filepath, scale=2):
    """Any global changes we want to make to figures before saving."""
    filepath = pathlib.Path(filepath)
    fig.write_image(filepath, scale=scale)
