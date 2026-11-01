"""Utilities for matching simulations with experiments."""

import logging
import os
import warnings
from copy import copy
from typing import IO

import numpy as np

from fairmd.lipids import FMDL_SIMU_PATH, progress
from fairmd.lipids.auxiliary import encode_canonical_yaml
from fairmd.lipids.core import System, SystemsCollection
from fairmd.lipids.experiment import Experiment, ExperimentCollection

logger = logging.getLogger("__name__")


class E2SComparator:
    """Compare experiment and simulation conditions."""

    @staticmethod
    def match_membrane_composition(sim_mf: dict, exp_mf: dict) -> bool:
        """Compare two membrane compositions given as molar fractions."""
        if not (set(sim_mf.keys()) & set(exp_mf.keys())):
            return False
        _sim_mf = copy(sim_mf)
        _exp_mf = copy(exp_mf)
        # Presence tolerance: remove molecules with very small fractions
        presence_tol = 0.01
        for key in set(sim_mf.keys()) - set(exp_mf.keys()):
            if sim_mf[key] < presence_tol:
                del _sim_mf[key]
        for key in set(exp_mf.keys()) - set(sim_mf.keys()):
            if exp_mf[key] < presence_tol:
                del _exp_mf[key]
        if set(_sim_mf.keys()) != set(_exp_mf.keys()):
            return False
        # SQRT-based tolerance
        sqrt_tolerance_molfraction = 0.03
        membrane_composition_ok = True
        for key in _sim_mf:
            if np.abs(np.sqrt(_exp_mf[key]) - np.sqrt(_sim_mf[key])) > sqrt_tolerance_molfraction:
                # more sensitive to small fraction changes and less to large fraction changes
                membrane_composition_ok = False
                break
        return membrane_composition_ok

    @staticmethod
    def match_solution_composition(sim_mf: dict, exp_mf: dict) -> bool:
        """Compare two solution compositions given as molar fractions."""
        if not sim_mf and not exp_mf:
            return True  # pure water
        if not (set(sim_mf.keys()) & set(exp_mf.keys())):
            return False
        _sim_mf = copy(sim_mf)
        _exp_mf = copy(exp_mf)
        presence_tol = 0.02
        for key in set(sim_mf.keys()) - set(exp_mf.keys()):
            if sim_mf[key] < presence_tol:
                del _sim_mf[key]
        for key in set(exp_mf.keys()) - set(sim_mf.keys()):
            if exp_mf[key] < presence_tol:
                del _exp_mf[key]
        if set(_sim_mf.keys()) != set(_exp_mf.keys()):
            return False
        solution_rel_log10_tol = 0.5
        return all(np.abs(np.log10(_exp_mf[key] / _sim_mf[key])) <= solution_rel_log10_tol for key in _sim_mf)

    @staticmethod
    def match_temperature(sim_t: float, exp_t: float) -> bool:
        """Check if two temperatures match within tolerance."""
        abs_tolerance_t = 2.5
        return np.abs(sim_t - exp_t) <= abs_tolerance_t

    @staticmethod
    def match_hydration(sim_h: float, exp_h: float) -> bool:
        """Check if two hydration levels match within tolerance."""
        if exp_h > 25 and sim_h > 25:
            return True
        lip_rel_hydration_tol = 0.15
        return np.abs(exp_h / sim_h - 1) <= lip_rel_hydration_tol


def _register_paired_experiments(simulation: System, experiments: list[Experiment]) -> None:
    """Add matched experiments to a simulation and sort their identifiers."""
    for experiment in experiments:
        if experiment.exptype == "OrderParameters":
            for lipid in experiment.data:
                if lipid not in simulation["EXPERIMENT"]["ORDERPARAMETER"]:
                    simulation["EXPERIMENT"]["ORDERPARAMETER"][lipid] = []
                simulation["EXPERIMENT"]["ORDERPARAMETER"][lipid].append(experiment.exp_id)
        elif experiment.exptype == "FormFactors":
            simulation["EXPERIMENT"]["FORMFACTOR"].append(experiment.exp_id)

    cur_exp = simulation["EXPERIMENT"]
    cur_exp["FORMFACTOR"].sort()
    for lipid in cur_exp["ORDERPARAMETER"]:
        cur_exp["ORDERPARAMETER"][lipid].sort()


def _find_pairs_for_simulation(experiments: ExperimentCollection, simulation: System) -> list[Experiment]:
    """Find experiment matches for one simulation and update its README data."""
    sim_lipids_mf = simulation.membrane_composition(basis="molar")
    try:
        sim_ions_mf = simulation.solution_composition(basis="molar")
        sim_hydr = simulation.get_hydration(basis="number")
    except ValueError:
        # implicit water - hydration is not supported currently => don't pair
        warnings.warn(
            f"Simulation {simulation['ID']} has implicit water; skipping experiment matching.",
            stacklevel=2,
        )
        return []

    paired_exps = []
    for experiment in experiments:
        exp_lipids_mf = experiment.membrane_composition(basis="molar")
        if not E2SComparator.match_membrane_composition(sim_lipids_mf, exp_lipids_mf):
            continue

        exp_ions_mf = experiment.solution_composition(basis="molar")
        if not E2SComparator.match_solution_composition(sim_ions_mf, exp_ions_mf):
            continue

        exp_hydr = experiment.get_hydration(basis="number")
        if not E2SComparator.match_hydration(sim_hydr, exp_hydr):
            continue

        if not E2SComparator.match_temperature(simulation["TEMPERATURE"], experiment["TEMPERATURE"]):
            continue

        paired_exps.append(experiment)

    _register_paired_experiments(simulation, paired_exps)

    return paired_exps


def pair_all_simulations(experiments: ExperimentCollection, simulations: SystemsCollection) -> list:
    """Find matching simulation-experiment pairs and update simulations' README data."""
    pairs = []
    for simulation in progress(simulations, desc="Simulation"):
        pairs.extend(pair_one_simulation(experiments, simulation))
    return pairs


def pair_one_simulation(experiments: ExperimentCollection, simulation: System) -> list[list[System | Experiment]]:
    """Pair one simulation with experiments and write its README to disk."""
    paired_exps = _find_pairs_for_simulation(experiments, simulation)
    outfile = os.path.join(FMDL_SIMU_PATH, simulation["path"], "README.yaml")
    with open(outfile, "w", encoding="utf-8") as f:
        f.write(encode_canonical_yaml(simulation.to_readme()))
    return [[simulation, experiment] for experiment in paired_exps]


def log_pairs(pairs, fd: IO[str]) -> None:
    """
    Write found correspondences into log file.

    pairs: [(Simulation, Experiment), ...]
    fd: file descriptor for writting into
    """
    for p in pairs:
        sim: System = p[0]
        exp: Experiment = p[1]

        sysn = sim["SYSTEM"]
        simp = sim["path"]

        expp = exp.path
        # The two DOI fields an experiment may give. The deprecated plain `DOI:`
        # key is deliberately not read here: renaming it is what the deprecation
        # asks for, and a log line is not a reason to keep it alive.
        expd = exp.readme.get("ARTICLE_DOI") or exp.readme.get("DATA_DOI") or "[no DOI]"

        fd.write(f"""
--------------------------------------------------------------------------------
Simulation:
 - {sysn}
 - {simp}
Experiment:
 - {expd}
 - {expp}""")
        # end for
    fd.write("""
--------------------------------------------------------------------------------
    \n""")
