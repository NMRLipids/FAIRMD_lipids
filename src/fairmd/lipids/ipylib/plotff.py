"""Plotting routines for form factors."""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure

import fairmd.lipids.analib.formfactor as ff
from fairmd.lipids.api import get_FF, get_quality
from fairmd.lipids.core import System
from fairmd.lipids.experiment import ExperimentCollection

from .style import fmdl_plot_style


def plot_simulation_FF(system: System) -> Figure:  # noqa: N802
    """Plot the simulated and experimental form factors for ``system``.

    NOTE: Currently, it plots only the first form factor experiment found in
    the system's metadata (if found).

    :return: The Matplotlib figure containing the form-factor plot.
    """
    print("DOI: ", system["DOI"])
    ff_quality = get_quality(system, experiment="FF")
    print("Form factor quality: ", ff_quality)

    ff_sim = np.asarray(get_FF(system), dtype=float)
    if len(system["EXPERIMENT"]["FORMFACTOR"]):
        form_factor = system["EXPERIMENT"]["FORMFACTOR"][0]
        ff_experiments = ExperimentCollection.load_from_data("FFExperiment")
        experiment = ff_experiments.get(form_factor)
        if experiment is None:
            msg = f"Form factor experiment {form_factor} not found in the database."
            raise FileNotFoundError(msg)
        ff_exp = np.asarray(experiment.data, dtype=float)
        scf = ff.calc_ff_scaling_distance(ff_exp, ff_sim)[0]
    else:
        ff_exp = None

    with plt.rc_context(fmdl_plot_style["common"]):
        fig, ax = plt.subplots()
        _plot_form_factor(ax, ff_sim, 1, "Simulation", "red")
        if ff_exp is not None:
            _plot_form_factor(ax, ff_exp, scf, "Experiment", "black")
        fig.tight_layout()
    return fig


def _plot_form_factor(ax: plt.Axes, ff_df: np.ndarray, scaling_factor: float, plot_label: str, plot_color: str) -> None:
    """:meta private:"""
    _df = ff_df.copy()
    _df[:, 1] *= scaling_factor
    ax.plot(_df[:, 0], _df[:, 1], label=plot_label, color=plot_color, linewidth=4.0)
    ax.set_xlabel(r"$q_{z} [Å^{-1}]$")
    ax.set_ylabel(r"$|F(q_{z})|$")
    ax.set_xlim([0, 0.69])
    ax.set_ylim([-10, 250])
    ax.legend(loc="upper right")
