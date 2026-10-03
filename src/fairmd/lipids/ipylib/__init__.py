"""Routines for IPython/Jupyter notebooks."""

from .plotff import plot_simulation_FF
from .plotop import plot_simulation_OP, plotSimulation

__all__ = ["plotSimulation", "plot_simulation_FF", "plot_simulation_OP"]
