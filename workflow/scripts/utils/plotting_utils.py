"""Configure Plotly templates and default settings."""

import plotly.graph_objects as go
import plotly.io as pio

# Preferred styles
pio.templates["scientific"] = go.layout.Template(
    layout=dict(
        margin=dict(l=30, r=30, t=30, b=30),
        autosize=True,
        xaxis=dict(showgrid=True),
        yaxis=dict(showgrid=True),
        title=dict(x=0.5, xanchor="center"),
    )
)
pio.templates.default = "simple_white+scientific"
