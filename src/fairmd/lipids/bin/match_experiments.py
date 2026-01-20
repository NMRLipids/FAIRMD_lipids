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

import yaml
from tqdm import tqdm

from fairmd.lipids import FMDL_SIMU_PATH
from fairmd.lipids.api import lipids_set
from fairmd.lipids.core import System, initialize_databank
from fairmd.lipids.experiment import Experiment, ExperimentCollection

logger = logging.getLogger("__name__")

# TODO: REMOVE IT COMPLETELY!!!
ions_list = ["POT", "SOD", "CLA", "CAL"]  # should contain names of all ions

LIP_CONC_REL_THRESHOLD = 0.15  # relative acceptable error for determination
# of the hydration in ssNMR


# TODO: derive from Simulation (if not to remove at all!)
class SearchSystem:
    system: dict
    idx_path: str

    def __init__(self, readme):
        self.system: System = readme
        self.idx_path = readme["path"]

    def get_lipids(self, molecules=lipids_set):
        """Return list of lipids"""
        lipids = [k for k in self.system["COMPOSITION"] if k in molecules]
        return lipids

    def get_ions(self, ions):
        """Return list of non-zero ions"""
        sim_ions = [k for k in self.system["COMPOSITION"] if k in ions]
        return sim_ions

    # fraction of each lipid with respect to total amount of lipids (only for lipids!)
    def molar_fraction(self, molecule, molecules=lipids_set) -> float:
        cmps = self.system["COMPOSITION"]
        number = sum(cmps[molecule]["COUNT"])
        all_counts = [i["COUNT"] for k, i in cmps.items() if k in molecules]
        return number / sum(map(sum, all_counts))

    # concentration of other molecules than lipids
    # change name to ionConcentration()

    def ion_conc(self, molecule, exp_counter_ions):
        lipids1 = self.get_lipids()
        c_water = 55.5
        n_water = self.system["COMPOSITION"]["SOL"]["COUNT"]
        try:
            n_molecule = self.system["COMPOSITION"][molecule]["COUNT"]  # number of ions
        except KeyError:
            n_molecule = 0

        lipids2 = []
        if exp_counter_ions and n_molecule != 0:
            for lipid in lipids1:
                if molecule in exp_counter_ions.keys() and lipid == exp_counter_ions[molecule]:
                    n_lipid = self.system["COMPOSITION"][lipid]["COUNT"]
                    lipids2.append(sum(n_lipid))

        n_molecule = n_molecule - sum(lipids2)
        c_molecule = (n_molecule * c_water) / n_water

        return c_molecule

    def total_lipid_conc(self):
        c_water = 55.5
        n_water = self.system["COMPOSITION"]["SOL"]["COUNT"]
        n_lipids = 0
        for lipid in self.get_lipids():
            try:
                n_lipids += sum(self.system["COMPOSITION"][lipid]["COUNT"])
            except KeyError as e:
                print(self.system)
                raise e
        try:
            if (n_water / n_lipids) > 25:
                tot_lipid_c = "full hydration"
            else:
                tot_lipid_c = (n_lipids * c_water) / n_water
        except ZeroDivisionError:
            logger.warning("Division by zero when determining lipid concentration!")
            print(self.system)
        return tot_lipid_c


##################


def load_simulations() -> list[SearchSystem]:
    """Generate the list of Simulation objects. Go through all README.yaml files."""
    systems = initialize_databank()
    simulations: list[SearchSystem] = []

    for system in systems:
        # conditions of exclusions
        try:
            if system["WARNINGS"]["NOWATER"]:
                continue
        except (KeyError, TypeError):
            pass

        simulations.append(SearchSystem(system))

    return simulations


def load_experiments(exp_type: str, all_experiments: ExperimentCollection) -> list[Experiment]:
    """Filter experiments from the collection by experiment type."""
    print(f"Filtering for {exp_type} experiments...")
    return [exp for exp in all_experiments if exp.exptype == exp_type]


def find_pairs_and_change_sims(experiments: list[Experiment], simulations: list[SearchSystem]):
    pairs = []
    for simulation in tqdm(simulations, desc="Simulation"):
        sim_lipids = simulation.get_lipids()
        sim_total_lipid_concentration = simulation.total_lipid_conc()
        sim_ions = simulation.get_ions(ions_list)
        t_sim = simulation.system["TEMPERATURE"]

        # calculate molar fractions from simulation
        sim_molar_fractions = {}
        for lipid in sim_lipids:
            sim_molar_fractions[lipid] = simulation.molar_fraction(lipid)

        for experiment in experiments:
            # check lipid composition matches the simulation
            exp_lipids = experiment.get_lipids()

            exp_total_lipid_concentration = experiment.readme["TOTAL_LIPID_CONCENTRATION"]
            exp_ions = experiment.get_ions(ions_list)
            exp_counter_ions = experiment.readme.get("COUNTER_IONS")

            # calculate simulation ion concentrations
            sim_concentrations = {}
            for molecule in ions_list:
                sim_concentrations[molecule] = simulation.ion_conc(molecule, exp_counter_ions)

            # continue if lipid compositions are the same
            if set(sim_lipids) == set(exp_lipids):
                # compare molar fractions
                mf_ok = 0
                for key in sim_lipids:
                    if (experiment.readme["MOLAR_FRACTIONS"][key] >= sim_molar_fractions[key] - 0.03) and (
                        experiment.readme["MOLAR_FRACTIONS"][key] <= sim_molar_fractions[key] + 0.03
                    ):
                        mf_ok += 1

                # compare ion concentrations
                c_ok = 0
                if set(sim_ions) == set(exp_ions):
                    for key in sim_ions:
                        if (experiment.readme["ION_CONCENTRATIONS"][key] >= sim_concentrations[key] - 0.05) and (
                            experiment.readme["ION_CONCENTRATIONS"][key] <= sim_concentrations[key] + 0.05
                        ):
                            c_ok += 1

                switch = 0

                if isinstance(exp_total_lipid_concentration, (int, float)) and isinstance(
                    sim_total_lipid_concentration,
                    (int, float),
                ):
                    if (
                        exp_total_lipid_concentration / sim_total_lipid_concentration > 1 - LIP_CONC_REL_THRESHOLD
                    ) and (exp_total_lipid_concentration / sim_total_lipid_concentration < 1 + LIP_CONC_REL_THRESHOLD):
                        switch = 1
                elif (
                    (type(exp_total_lipid_concentration) is str)
                    and (type(sim_total_lipid_concentration) is str)
                    and (exp_total_lipid_concentration == sim_total_lipid_concentration)
                ):
                    switch = 1

                if switch:
                    # check temperature +/- 2 degrees
                    t_exp = experiment.readme["TEMPERATURE"]

                    if (
                        (mf_ok == len(sim_lipids))
                        and (c_ok == len(sim_ions))
                        and (t_exp > float(t_sim) - 2.5)
                        and (t_exp < float(t_sim) + 2.5)
                    ):
                        # !we found the match!
                        pairs.append([simulation, experiment])

                        # Add path to experiment into simulation README.yaml
                        # many experiment entries can match to same simulation
                        if experiment.exptype == "OrderParameters":
                            for lipid in experiment.data:
                                if lipid not in simulation.system["EXPERIMENT"]["ORDERPARAMETER"]:
                                    simulation.system["EXPERIMENT"]["ORDERPARAMETER"][lipid] = []
                                simulation.system["EXPERIMENT"]["ORDERPARAMETER"][lipid].append(experiment.exp_id)
                        elif experiment.exptype == "FormFactors":
                            simulation.system["EXPERIMENT"]["FORMFACTOR"].append(experiment.exp_id)
                    else:
                        continue

        # sorting experiment lists to keep experimental order strict
        cur_exp = simulation.system["EXPERIMENT"]
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
        sim: SearchSystem = p[0]
        exp: Experiment = p[1]

        sysn = sim.system["SYSTEM"]
        simp = sim.idx_path

        expp = exp.dataPath
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
    simulations = load_simulations()

    # clear all EXPERIMENT sections in all simulations
    for simulation in simulations:
        simulation.system["EXPERIMENT"] = {}
        simulation.system["EXPERIMENT"]["ORDERPARAMETER"] = {}
        simulation.system["EXPERIMENT"]["FORMFACTOR"] = []
        for lipid in simulation.get_lipids():
            simulation.system["EXPERIMENT"]["ORDERPARAMETER"][lipid] = []

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
        outfile_dict = os.path.join(FMDL_SIMU_PATH, simulation.idx_path, "README.yaml")
        with open(outfile_dict, "w") as f:
            if "path" in simulation.system:
                del simulation.system["path"]
            yaml.dump(simulation.system.readme, f, sort_keys=False, allow_unicode=True)

    print("Found order parameter data for " + str(len(pairs_op)) + " pairs")
    print("Found form factor data for " + str(len(pairs_ff)) + " pairs")


if __name__ == "__main__":
    match_experiments()
