#!/usr/bin/env python3
"""
Perform comparison of experiments and simulations.

The script compares according to **EXPERIMENT** field inside :ref:`the simulation README.yaml file <readmesimu>`.
In the standard protocol, it should be run *after* :ref:`fmdl_match_experiments <match_experiments_py>`.

**Usage:**

.. code-block:: console

    fmdl_evaluate_quality

No arguments are needed.
"""

import json
import os

import numpy as np
import yaml

import fairmd.lipids.quality as qq
from fairmd.lipids import FMDL_EXP_PATH, FMDL_SIMU_PATH
from fairmd.lipids.jsonEncoders import CompactJSONEncoder


def evaluate_quality():
    simulations = qq.loadSimulations()

    EvaluatedOPs = 0
    EvaluatedFFs = 0

    for simulation in simulations:
        # save OP quality and FF quality here
        DATAdir = os.path.join(FMDL_SIMU_PATH, simulation.idx_path)
        print("Analyzing: ", DATAdir)

        # Order Parameters
        system_quality = {}
        for lipid1 in simulation.get_lipids():
            print(f"\nEvaluating order parameter quality of simulation data in {simulation.idx_path}")

            OP_data_lipid = {}
            # convert elements to float because in some files the elements are strings
            try:
                for key, _ in simulation.op_data[lipid1].items():
                    OP_array = [float(x) for x in simulation.op_data[lipid1][key][0]]
                    OP_data_lipid[key] = OP_array
            except Exception:
                continue

            # go through file paths in simulation.readme['EXPERIMENT']
            fragment_qual_dict = {}
            data_dict = {}

            for doi, path in simulation.system["EXPERIMENT"]["ORDERPARAMETER"][lipid1].items():
                print(
                    f"Evaluating {lipid1} lipid using experimental data from"
                    f"{doi} in {FMDL_EXP_PATH}/OrderParameters/{path}",
                )

                print(doi)
                OP_qual_data = {}
                # get readme file of the experiment
                exp_fpath = os.path.join(FMDL_EXP_PATH, "OrderParameters", path)
                print("Experimental data available at " + exp_fpath)

                READMEfilepathExperiment = os.path.join(exp_fpath, "README.yaml")
                experiment = qq.Experiment()
                with open(READMEfilepathExperiment) as yaml_file_exp:
                    readme_exp = yaml.load(yaml_file_exp, Loader=yaml.FullLoader)
                    experiment.readme = readme_exp

                exp_op_fpath = os.path.join(exp_fpath, lipid1 + "_Order_Parameters.json")
                exp_op_data = {}
                try:
                    with open(exp_op_fpath) as json_file:
                        exp_op_data = json.load(json_file)
                except FileNotFoundError:
                    print(f"Experimental order parameter data do not exist for lipid {lipid1}.")
                    continue
                except Exception as e:
                    raise RuntimeError(f"Unexpected error during loading {exp_op_fpath}") from e

                exp_error = 0.02

                for key in OP_data_lipid:
                    OP_array = OP_data_lipid[key].copy()
                    try:
                        OP_exp = exp_op_data[key][0][0]
                    except KeyError:
                        continue
                    else:
                        if not np.isnan(OP_exp):
                            OP_sim = OP_array[0]
                            op_sim_STEM = OP_array[2]
                            # changing to use shitness(TM) scale.
                            # This code needs to be cleaned
                            op_quality = qq.prob_S_in_g(OP_exp, exp_error, OP_sim, op_sim_STEM)
                            OP_array.append(OP_exp)
                            OP_array.append(exp_error)  # hardcoded!!! 0.02 for all exps
                            OP_array.append(op_quality)
                    OP_qual_data[key] = OP_array

                # save qualities of simulation-vs-experiment into a dictionary
                data_dict[doi] = OP_qual_data

                # calculate quality for molecule fragments headgroup, sn-1, sn-2
                fragments = qq.get_fragments(simulation.system.content[lipid1].mapping_dict)
                fragment_qual_dict[doi] = qq.fragmentQuality(fragments, exp_op_data, OP_data_lipid)

            try:
                fragment_quality_output = qq.fragmentQualityAvg(lipid1, fragment_qual_dict, fragments)
            except Exception:
                print("no fragment quality")
                fragment_quality_output = {}

            try:
                system_quality[lipid1] = fragment_quality_output
            except Exception:
                print("no system quality")
                system_quality[lipid1] = {}

            fragment_quality_file = os.path.join(DATAdir, lipid1 + "_FragmentQuality.json")

            FGout = False
            for FG in fragment_quality_output:
                # print(FG,fragment_quality_output[FG])
                if np.isnan(fragment_quality_output[FG]):
                    continue
                if fragment_quality_output[FG] > 0:
                    FGout = True
            if FGout:
                # write fragment qualities into a file for a molecule
                with open(fragment_quality_file, "w") as f:
                    json.dump(fragment_quality_output, f)

            # write into the OrderParameters_quality.json quality data file
            outfile1 = os.path.join(DATAdir, lipid1 + "_OrderParameters_quality.json")
            try:
                with open(outfile1, "w") as f:
                    json.dump(data_dict, f, cls=CompactJSONEncoder)
            except Exception:
                pass

        system_qual_output = qq.systemQuality(system_quality, simulation)
        # make system quality file

        outfile2 = os.path.join(DATAdir, "SYSTEM_quality.json")
        SQout = False
        for SQ in system_qual_output:
            if system_qual_output[SQ] > 0:
                SQout = True
        if SQout:
            with open(outfile2, "w") as f:
                json.dump(system_qual_output, f)
            print("Order parameter quality evaluated for " + simulation.idx_path)
            EvaluatedOPs += 1
            print()

        ###############################################################################
        # Form factor quality

        expFFpath = simulation.system["EXPERIMENT"]["FORMFACTOR"]
        expFFdata = {}
        if len(expFFpath) > 0:
            expFFpath_full = os.path.join(FMDL_EXP_PATH, "FormFactors", expFFpath)
            for _, _, files in os.walk(expFFpath_full):
                for filename in files:
                    filepath = os.path.join(expFFpath_full, filename)
                    if filename.endswith(".json"):
                        with open(filepath) as json_file:
                            expFFdata = json.load(json_file)

        simFFdata = simulation.ff_data

        if len(expFFpath) > 0 and len(simFFdata) > 0:
            ffQuality = qq.formfactorQuality(simFFdata, expFFdata)
            outfile3 = os.path.join(DATAdir, "FormFactorQuality.json")
            with open(outfile3, "w") as f:
                json.dump(ffQuality, f)
            EvaluatedFFs += 1
            print("Form factor quality evaluated for ", DATAdir)
        else:
            ffQuality = 0

    print("The number of systems with evaluated order parameters:", EvaluatedOPs)
    print("The number of systems with evaluated form factors:", EvaluatedFFs)


if __name__ == "__main__":
    evaluate_quality()
