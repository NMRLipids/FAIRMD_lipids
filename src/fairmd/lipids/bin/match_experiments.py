#!/usr/bin/env python3
"""Command-line entry point for matching simulations and experiments."""

from fairmd.lipids.core import initialize_databank
from fairmd.lipids.experiment import ExperimentCollection
from fairmd.lipids.matcher import log_pairs, pair_all_simulations


def match_experiments() -> None:
    """Run matching for all databank simulations and experiment types."""
    simulations = initialize_databank()

    # clear all EXPERIMENT sections in all simulations
    for simulation in simulations:
        simulation["EXPERIMENT"] = {}
        simulation["EXPERIMENT"]["ORDERPARAMETER"] = {}
        simulation["EXPERIMENT"]["FORMFACTOR"] = []
        for lipid in simulation.lipids:
            simulation["EXPERIMENT"]["ORDERPARAMETER"][lipid] = []

    with open("search-databank-pairs.log", "w") as logf:
        print("Scanning simulation-experiment pairs among order parameter experiments.")
        exps = ExperimentCollection.load_from_data("OPExperiment")
        print(f"{len(exps)} OP experiments loaded.")
        pairs_op = pair_all_simulations(exps, simulations)
        logf.write("=== OP PAIRS ===\n")
        log_pairs(pairs_op, logf)

        exps = ExperimentCollection.load_from_data("FFExperiment")
        print(f"{len(exps)} FF experiments loaded.")
        print("Scanning simulation-experiment pairs among form factor experiments.")
        pairs_ff = pair_all_simulations(exps, simulations)
        logf.write("=== FF PAIRS ===\n")
        log_pairs(pairs_ff, logf)

    print("Found order parameter data for " + str(len(pairs_op)) + " pairs")
    print("Found form factor data for " + str(len(pairs_ff)) + " pairs")


if __name__ == "__main__":
    match_experiments()
