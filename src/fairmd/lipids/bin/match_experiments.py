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
"""

import logging
import os
from typing import IO

import yaml
from tqdm import tqdm

from fairmd.lipids import FMDL_EXP_PATH, FMDL_SIMU_PATH
from fairmd.lipids.core import System, initialize_databank
from fairmd.lipids.databankLibrary import lipids_set

logger = logging.getLogger("__name__")

# TODO: move ions list into Data
ions_list = ["POT", "SOD", "CLA", "CAL"]  # should contain names of all ions

LIP_CONC_REL_THRESHOLD = 0.15  # relative acceptable error for determination
# of the hydration in ssNMR


class SearchSystem:
    system: dict
    idx_path: str

    def __init__(self, readme):
        self.system: System = readme
        self.idx_path = readme["path"]

    def get_lipids(self, molecules=lipids_set):
        lipids = [k for k in self.system["COMPOSITION"] if k in molecules]
        return lipids

    def get_ions(self, ions):
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


class Experiment:
    def __init__(self, readme: dict, molname: str, data_path: str, exptype: str):
        self.readme = readme
        self.molname = molname  # <- the dictionary about existence of data files
        self.dataPath = data_path  # .... for particular lipids
        self.exptype = exptype

    def get_lipids(self, molecules=lipids_set) -> list[str]:
        lipids = [k for k in self.readme["MOLAR_FRACTIONS"] if k in molecules]
        return lipids

    def get_ions(self, ions) -> list[str]:
        exp_ions: list[str] = []

        for key in ions:
            try:
                if self.readme["ION_CONCENTRATIONS"][key] != 0:
                    exp_ions.append(key)
            except KeyError:
                continue
            try:
                if key in self.readme["COUNTER_IONS"]:
                    exp_ions.append(key)
            except (TypeError, KeyError):
                continue
        return exp_ions


def load_simulations() -> list[SearchSystem]:
    """
    Generates the list of Simulation objects. Go through all README.yaml files.
    """
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


def load_experiments(exp_type: str) -> list[Experiment]:
    """
    Loops over the experiment entries in the experiment databank and read experiment
    readme and order parameter files into objects.
    """
    if exp_type == "OrderParameters":
        data_file = "_Order_Parameters.json"
    elif exp_type == "FormFactors":
        data_file = "_FormFactor.json"
    else:
        raise NotImplementedError("Only OrderParameters and FormFactors types are implemented.")

    print("Build experiments [%s] index..." % exp_type, end="")
    rm_idx = []

    path = os.path.join(FMDL_EXP_PATH, exp_type)
    for subdir, _, files in os.walk(path):
        for fn in files:
            if fn == "README.yaml":
                rm_idx.append(subdir)
    print("%d READMEs loaded." % len(rm_idx))

    print("Loading data for each experiment.")
    experiments: list[Experiment] = []
    for subdir in tqdm(rm_idx, desc="Experiment"):
        try:
            exp_readme_fp = os.path.join(subdir, "README.yaml")
            with open(exp_readme_fp) as yaml_file_exp:
                exp_readme = yaml.load(yaml_file_exp, Loader=yaml.FullLoader)
        except (FileNotFoundError, PermissionError):
            logger.warning(f"Problems while accessing README.yaml in: {subdir}")
            continue

        for fname in os.listdir(subdir):
            if fname.endswith(data_file):
                molecule_name = ""
                if exp_type == "OrderParameters":
                    molecule_name = fname.replace(data_file, "")
                elif exp_type == "FormFactors":
                    molecule_name = "system"
                experiments.append(Experiment(exp_readme, molecule_name, subdir, exp_type))

    return experiments


def find_pairs(experiments: list[Experiment], simulations: list[SearchSystem]):
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
            exp_counter_ions = experiment.readme["COUNTER_IONS"]

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
                elif (type(exp_total_lipid_concentration) is str) and (type(sim_total_lipid_concentration) is str):
                    if exp_total_lipid_concentration == sim_total_lipid_concentration:
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
                        exp_doi = experiment.readme["DOI"]
                        exp_path = os.path.relpath(
                            experiment.dataPath,
                            start=os.path.join(FMDL_EXP_PATH, experiment.exptype),
                        )
                        if experiment.exptype == "OrderParameters":
                            lipid = experiment.molname
                            simulation.system["EXPERIMENT"]["ORDERPARAMETER"][lipid][exp_doi] = exp_path
                        elif experiment.exptype == "FormFactors":
                            simulation.system["EXPERIMENT"]["FORMFACTOR"] = exp_path
                    else:
                        continue

        # sorting experiment lists to keep experimental order strict
        cur_exp = simulation.system["EXPERIMENT"]
        for _lipid in cur_exp["ORDERPARAMETER"]:
            unsort_dict = cur_exp["ORDERPARAMETER"][_lipid].copy()
            if not len(unsort_dict):
                continue
            sort_dict = dict(sorted(unsort_dict.items()))
            cur_exp["ORDERPARAMETER"][_lipid] = sort_dict.copy()

        outfile_dict = os.path.join(FMDL_SIMU_PATH, simulation.idx_path, "README.yaml")
        with open(outfile_dict, "w") as f:
            if "path" in simulation.system.keys():
                del simulation.system["path"]
            yaml.dump(simulation.system.readme, f, sort_keys=False, allow_unicode=True)

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
        expd = exp.readme["DOI"]

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


def match_experiments():
    """
    Main program function. Not for exporting.
    """
    simulations = load_simulations()

    # clear all EXPERIMENT sections in all simulations
    # TODO: check if EXPERIMENT section changed and trigger the action!
    for simulation in simulations:
        simulation.system["EXPERIMENT"] = {}
        simulation.system["EXPERIMENT"]["ORDERPARAMETER"] = {}
        simulation.system["EXPERIMENT"]["FORMFACTOR"] = {}
        for lipid in simulation.get_lipids():
            simulation.system["EXPERIMENT"]["ORDERPARAMETER"][lipid] = {}

        readme_path = os.path.join(FMDL_SIMU_PATH, simulation.idx_path, "README.yaml")
        with open(readme_path, "w") as f:
            yaml.dump(simulation.system.readme, f, sort_keys=False, allow_unicode=True)

    experiments_op = load_experiments("OrderParameters")
    experiments_ff = load_experiments("FormFactors")

    # Pair each simulation with an experiment with the closest matching temperature
    # and composition
    with open("search-databank-pairs.log", "w") as logf:
        print("Scanning simulation-experiment pairs among order parameter experiments.")
        pairs_op = find_pairs(experiments_op, simulations)
        logf.write("=== OP PAIRS ===\n")
        log_pairs(pairs_op, logf)
        print("Scanning simulation-experiment pairs among form factor experiments.")
        pairs_ff = find_pairs(experiments_ff, simulations)
        logf.write("=== FF PAIRS ===\n")
        log_pairs(pairs_ff, logf)

    """
    for pair in pairsFF:
        print('#################')
        print(pair[0].readme)
        print(pair[0].indexingPath)
        print("#")
        print(pair[1].readme)
    """

    print("Found order parameter data for " + str(len(pairs_op)) + " pairs")
    print("Found form factor data for " + str(len(pairs_ff)) + " pairs")


if __name__ == "__main__":
    match_experiments()
