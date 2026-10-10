"""Plotting style configurations"""

# Define default Matplotlib rc parameters shared by the plotting routines.
fmdl_plot_style = {
    # Every key in this mapping must be a valid Matplotlib rcParam because it is
    # passed directly to ``plt.rc_context``.
    "common": {
        "font.size": 13,
        "figure.figsize": (8.5, 5.5),
        "figure.dpi": 120,
        "axes.labelsize": 16,
        "xtick.labelsize": 13,
        "ytick.labelsize": 13,
        "xtick.major.width": 1.2,
        "ytick.major.width": 1.2,
        "xtick.major.size": 6,
        "ytick.major.size": 6,
    },
    # OP for simulation
    "simulation": {
        "fmt": "s",
        "label": "Simulation",
        "color": "red",
        "markersize": 9,
        "markeredgecolor": "black",
        "markeredgewidth": 0.7,
        "capsize": 2,
    },
    # OP for experiment
    "experimental": {
        "fmt": "o",
        "label": "Experimental",
        "color": "black",
        "markersize": 9,
        "markeredgecolor": "gray",
        "markeredgewidth": 0.7,
        "capsize": 2,
    },
}
