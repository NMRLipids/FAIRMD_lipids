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
    the system's metadata.

    :return: The Matplotlib figure containing the form-factor plot.
    """
    print("DOI: ", system["DOI"])
    ff_quality = get_quality(system, experiment="FF")
    print("Form factor quality: ", ff_quality)

    ff_experiments = ExperimentCollection.load_from_data("FFExperiment")
    ff_exp = None
    for form_factor in system["EXPERIMENT"]["FORMFACTOR"]:
        experiment = ff_experiments.get(form_factor)
        if experiment is not None:
            ff_exp = experiment.data
            break
    if ff_exp is None:
        msg = "No form factor experiment was found"
        raise FileNotFoundError(msg)
    ff_exp = np.asarray(ff_exp, dtype=float)
    ff_sim = np.asarray(get_FF(system), dtype=float)
    scf = ff.calc_ff_scaling_distance(ff_exp, ff_sim)[0]

    plt.rcParams.update({"font.size": fmdl_plot_style["font.size"]})
    figure, _ = plt.subplots(
        figsize=fmdl_plot_style["figure.figsize"],
        dpi=fmdl_plot_style["figure.dpi"],
    )
    plotFormFactor(ff_sim, 1, "Simulation", "red")
    plotFormFactor(ff_exp, scf, "Experiment", "black")
    return figure


def plotFormFactor(exp_form_factor, k, legend, plot_color):  # noqa: N802
    """:meta private:"""
    x_vals = []
    y_vals = []
    for i in exp_form_factor:
        x_vals.append(i[0])
        y_vals.append(k * i[1])
    plt.plot(x_vals, y_vals, label=legend, color=plot_color, linewidth=4.0)
    plt.xlabel(r"$q_{z} [Å^{-1}]$", size=20)
    plt.ylabel(r"$|F(q_{z})|$", size=20)
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.xlim([0, 0.69])
    plt.ylim([-10, 250])
    plt.legend(loc="upper right")
    plt.tight_layout()
