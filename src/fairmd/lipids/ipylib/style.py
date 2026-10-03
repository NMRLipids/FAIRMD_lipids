"""Plotting style configurations"""

# Define a default plotting style for OP
fmdl_plot_style = {
    "font.size": 13,
    "figure.figsize": (8.5, 5.5),
    "figure.dpi": 120,
    "label_size": 16,
    "tick_size": 13,
    "tick_width": 1.2,
    "tick_length": 6,
    "simulation": {
        "fmt": "s",
        "label": "Simulation",
        "color": "red",
        "markersize": 9,
        "markeredgecolor": "black",
        "markeredgewidth": 0.7,
        "capsize": 2,
    },
    "experimental": {
        "fmt": "o",
        "label": "Experimental",
        "color": "blue",
        "markersize": 9,
        "markeredgecolor": "black",
        "markeredgewidth": 0.7,
        "capsize": 2,
    },
}
