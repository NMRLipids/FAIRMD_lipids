"""Implementation of plotting routines for IPython/Jupyter notebooks."""

import warnings

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.figure import Figure

from fairmd.lipids.api import get_OP
from fairmd.lipids.auxiliary.opconvertor import build_nice_OPdict
from fairmd.lipids.core import System
from fairmd.lipids.experiment import ExperimentCollection
from fairmd.lipids.molecules import Lipid

from .plotff import plot_simulation_FF
from .style import fmdl_plot_style


def _plot_general_lipid(op_nice_sim: dict, op_nice_exp: dict | None, lipid_obj: Lipid) -> tuple[Figure, Figure, Figure]:
    """Plot simulation and experimental order parameters by fragment."""

    def _prep_df(nicedic: dict) -> dict:
        for frag, data in nicedic.items():
            _df = pd.DataFrame(data)
            if _df.empty:
                nicedic[frag] = None
                continue
            _df.index = _df.index.astype(int)
            _df = _df.sort_index()
            _df["STD"] = pd.to_numeric(_df["STD"], errors="coerce")  # nan-ificate
            nicedic[frag] = _df
        return nicedic

    sim_df_dict = _prep_df(op_nice_sim)
    exp_df_dict = _prep_df(op_nice_exp) if op_nice_exp is not None else {}

    plt.rcParams.update({"font.size": fmdl_plot_style["font.size"]})
    fig_list = []
    for frag in sim_df_dict:
        simdf = sim_df_dict[frag]
        expdf = exp_df_dict.get(frag, None)
        figure, axis = plt.subplots(
            figsize=fmdl_plot_style["figure.figsize"],
            dpi=fmdl_plot_style["figure.dpi"],
        )
        axis.set_title(f"{lipid_obj.name} : {frag}", fontsize=fmdl_plot_style["label_size"])
        axis.errorbar(
            simdf["C"],
            simdf["OP"],
            yerr=simdf["STD"],
            **fmdl_plot_style["simulation"],
        )
        if expdf is not None:
            axis.errorbar(
                expdf["C"],
                expdf["OP"],
                yerr=expdf["STD"],
                **fmdl_plot_style["experimental"],
            )
        axis.set_xticks(simdf.C)
        axis.set_xlabel("Carbon", fontsize=fmdl_plot_style["label_size"])
        axis.set_ylabel(r"$S_{CH}$", fontsize=fmdl_plot_style["label_size"])
        axis.tick_params(
            axis="both",
            which="major",
            labelsize=fmdl_plot_style["tick_size"],
            width=fmdl_plot_style["tick_width"],
            length=fmdl_plot_style["tick_length"],
        )
        figure.tight_layout()
        fig_list.append(figure)
    return tuple(fig_list)


def plot_regular_phospholipid_op(op_sim: dict, op_exp: dict | None, lipid_obj: Lipid) -> tuple[Figure, Figure, Figure]:
    """Plot simulation and experimental order parameters by fragment."""

    def _unite_headgroup(op_dict: dict) -> dict:
        """Unite glycerol + head into head."""
        nicedic = build_nice_OPdict(op_dict, lipid_obj)
        nicedic["head"] = nicedic.get("glycerol backbone", []) + nicedic.get("headgroup", [])
        nicedic.pop("glycerol backbone", None)
        nicedic.pop("headgroup", None)
        return nicedic

    op_sim_nice = _unite_headgroup(op_sim)
    op_exp_nice = _unite_headgroup(op_exp) if op_exp else None

    return _plot_general_lipid(op_sim_nice, op_exp_nice, lipid_obj)


def plot_simulation_OP(system: System, lipid: str) -> tuple[Figure, Figure, Figure]:  # noqa: N802
    """Plot simulated and experimental C-H bond order parameters."""
    op_sim = get_OP(system).get(lipid)
    if op_sim is None:
        msg = f"Order parameter data not found for {lipid}"
        raise FileNotFoundError(msg)

    op_exp = {}
    op_experiments = ExperimentCollection.load_from_data("OPExperiment")
    for exp_op_id in system["EXPERIMENT"]["ORDERPARAMETER"][lipid]:
        experiment = op_experiments.get(exp_op_id)
        if experiment is not None:
            op_exp.update(experiment.data.get(lipid, {}))

    return plot_regular_phospholipid_op(op_sim, op_exp, system.lipids[lipid])


def plotSimulation(system: System, lipid: str) -> None:  # noqa: N802
    """Plot form factors and order parameters (deprecated)."""
    warnings.warn(
        "plotSimulation is deprecated; use plot_simulation_FF and plot_simulation_OP instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    plot_simulation_FF(system)
    plot_simulation_OP(system, lipid)
