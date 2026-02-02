#!/usr/bin/env python3
r"""
Match simulations with experiments in the databank.

Script goes through all simulations and experiments in the databank and finds
pairs of simulations and experiments that match in composition, temperature and
other conditions. The found pairs are written into the simulation :ref:`README.yaml <readmesimu>`
files and into a log file.

**Usage:**

.. code-block:: console

    fmdl_match_experiments

No arguments are needed.
TODO: check if EXPERIMENT section changed and trigger the action!
"""

import logging
import os
from typing import IO

import numpy as np
import yaml
from tqdm import tqdm

from fairmd.lipids import FMDL_SIMU_PATH
from fairmd.lipids.core import System, SystemsCollection, initialize_databank
from fairmd.lipids.experiment import Experiment, ExperimentCollection

logger = logging.getLogger("__name__")

# Helper functions for matching criteria
# TODO: remove to a module and TEST!


def _match_membrane_composition(sim_mf: dict, exp_mf: dict) -> bool:
    """Compare two membrane compositions given as molar fractions."""
    if set(sim_mf.keys()) != set(exp_mf.keys()):
        return False
    abs_tolerance_molfraction = 0.03
    # TODO: BAD! use relative threshold instead!
    membrane_composition_ok = True
    for key in sim_mf:
        if np.abs(exp_mf[key] - sim_mf[key]) > abs_tolerance_molfraction:
            membrane_composition_ok = False
            break
    return membrane_composition_ok


def _match_solution_composition(sim_mf: dict, exp_mf: dict) -> bool:
    """Compare two solution compositions given as molar fractions."""
    if set(sim_mf.keys()) != set(exp_mf.keys()):
        return False
    # BAD! use relative threshold instead!
    # BAD! use logarithmic scale instead!
    solution_rel_log10_tol = 0.5
    solution_composition_ok = True
    for key in sim_mf:
        if np.abs(np.log10(exp_mf[key] / sim_mf[key])) > solution_rel_log10_tol:
            solution_composition_ok = False
            break
    return solution_composition_ok


def _match_temperature(sim_t: float, exp_t: float) -> bool:
    """Check if two temperatures match within tolerance."""
    abs_tolerance_t = 2.5
    return np.abs(sim_t - exp_t) <= abs_tolerance_t


def _match_hydration(sim_h: float, exp_h: float) -> bool:
    """Check if two hydration levels match within tolerance."""
    if exp_h > 25 and sim_h > 25:  # both not full hydration
        return True
    lip_rel_hydration_tol = 0.15  # relative acceptable error for determination
    return np.abs(exp_h / sim_h - 1) <= lip_rel_hydration_tol


def find_pairs_and_change_sims(experiments: ExperimentCollection, simulations: SystemsCollection):
    """Find matching simulation-experiment pairs and update simulations' README data."""
    pairs = []
    for simulation in tqdm(simulations, desc="Simulation"):
        sim_lipids_mf = simulation.membrane_composition(basis="molar")
        try:
            sim_ions_mf = simulation.solution_composition(basis="molar")
            sim_hydr = simulation.get_hydration(basis="number")
        except ValueError:
            # implicit water - hydration is not supported currently => don't pair
            continue

        for experiment in experiments:
            exp_lipids_mf = experiment.membrane_composition(basis="molar")
            if not _match_membrane_composition(sim_lipids_mf, exp_lipids_mf):
                continue

            exp_ions_mf = experiment.solution_composition(basis="molar")
            if not _match_solution_composition(sim_ions_mf, exp_ions_mf):
                continue

            exp_hydr = experiment.get_hydration(basis="number")
            if not _match_hydration(sim_hydr, exp_hydr):
                continue

            if not _match_temperature(simulation["TEMPERATURE"], experiment["TEMPERATURE"]):
                continue

            pairs.append([simulation, experiment])

            # Add path to experiment into simulation README.yaml
            # many experiment entries can match to same simulation
            if experiment.exptype == "OrderParameters":
                for lipid in experiment.data:
                    if lipid not in simulation["EXPERIMENT"]["ORDERPARAMETER"]:
                        simulation["EXPERIMENT"]["ORDERPARAMETER"][lipid] = []
                    simulation["EXPERIMENT"]["ORDERPARAMETER"][lipid].append(experiment.exp_id)
            elif experiment.exptype == "FormFactors":
                simulation["EXPERIMENT"]["FORMFACTOR"].append(experiment.exp_id)

        # sorting experiment lists to keep experimental order strict
        cur_exp = simulation["EXPERIMENT"]
        cur_exp["FORMFACTOR"].sort()
        for _lipid in cur_exp["ORDERPARAMETER"]:
            cur_exp["ORDERPARAMETER"][_lipid].sort()

    return pairs


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
        expd = exp.readme.get("ARTICLE_DOI", "[no article DOI]")

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


def match_experiments() -> None:
    """Do main program work. Not for exporting."""
    simulations = initialize_databank()

    # clear all EXPERIMENT sections in all simulations
    for simulation in simulations:
        simulation["EXPERIMENT"] = {}
        simulation["EXPERIMENT"]["ORDERPARAMETER"] = {}
        simulation["EXPERIMENT"]["FORMFACTOR"] = []
        for lipid in simulation.lipids:
            simulation["EXPERIMENT"]["ORDERPARAMETER"][lipid] = []

    # Pair each simulation with an experiment with the closest matching temperature
    # and composition
    with open("search-databank-pairs.log", "w") as logf:
        print("Scanning simulation-experiment pairs among order parameter experiments.")
        exps = ExperimentCollection.load_from_data("OPExperiment")
        print(f"{len(exps)} OP experiments loaded.")
        pairs_op = find_pairs_and_change_sims(exps, simulations)
        logf.write("=== OP PAIRS ===\n")
        log_pairs(pairs_op, logf)

        exps = ExperimentCollection.load_from_data("FFExperiment")
        print(f"{len(exps)} FF experiments loaded.")
        print("Scanning simulation-experiment pairs among form factor experiments.")
        pairs_ff = find_pairs_and_change_sims(exps, simulations)
        logf.write("=== FF PAIRS ===\n")
        log_pairs(pairs_ff, logf)

    # save changed simulations
    for simulation in tqdm(simulations, "Saving READMEs"):
        outfile_dict = os.path.join(FMDL_SIMU_PATH, simulation["path"], "README.yaml")
        with open(outfile_dict, "w") as f:
            if "path" in simulation:
                del simulation["path"]
            yaml.dump(simulation.readme, f, sort_keys=False, allow_unicode=True)

    print("Found order parameter data for " + str(len(pairs_op)) + " pairs")
    print("Found form factor data for " + str(len(pairs_ff)) + " pairs")


if __name__ == "__main__":
    match_experiments()
