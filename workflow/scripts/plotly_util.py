"""Configure Plotly templates and default settings."""

import datetime
import pathlib

import plotly.graph_objects as go
import plotly.io as pio

pio.templates["pds"] = go.layout.Template(
    layout=dict(
        margin=dict(l=30, r=30, t=30, b=30),
        autosize=True,
        width=600,
        height=400,
        xaxis=dict(showgrid=True),
        yaxis=dict(showgrid=True),
        title=dict(x=0.5, xanchor="center"),
    )
)
pio.templates.default = "simple_white+pds"


def save_plotly_fig(fig, filepath, scale=2):
    """Any global changes we want to make to figures before saving."""
    filepath = pathlib.Path(filepath)
    fig.write_image(filepath, scale=scale)
