#!/usr/bin/env python3
"""
Creates different types of ranking lists inside the Databank.

It ranks simulations based on their quality against experiments. The ranking lists are stored in
``{FMDL_DATA_PATH}/Ranking/`` folder in JSON format.
The lists can be shown with the ``fairmd.lipids.plottings`` module.

**Usage:**

.. code-block:: console

    fmdl_make_ranking

No arguments are needed.
"""

import json
import os

from fairmd.lipids import FMDL_DATA_PATH, FMDL_SIMU_PATH
from fairmd.lipids.core import initialize_databank
from fairmd.lipids.databankLibrary import lipids_set
from fairmd.lipids.jsonEncoders import CompactJSONEncoder


def make_ranking():
    systems = initialize_databank()

    # ---- Making list of qualities
    qualities = []
    for system in systems:
        path = os.path.join(FMDL_SIMU_PATH, system["path"])

        total_quality_file_path = os.path.join(path, "SYSTEM_quality.json")

        fragment_Q = {}  # noqa: N816
        for lipid in system["COMPOSITION"]:
            quality_file = os.path.join(path, lipid + "_FragmentQuality.json")
            try:
                with open(quality_file) as json_file:
                    fragment_Q[lipid] = json.load(json_file)
            except Exception:
                continue

        if os.path.isfile(total_quality_file_path):
            with open(total_quality_file_path) as json_file:
                fragment_Q["TotalQuality"] = json.load(json_file)

        ff_quality_fpath = os.path.join(path, "FormFactorQuality.json")
        if os.path.isfile(ff_quality_fpath):
            with open(ff_quality_fpath) as json_file:
                FFq = json.load(json_file)
            try:
                fragment_Q["TotalQuality"]
            except KeyError:
                fragment_Q["TotalQuality"] = {}
            try:
                fragment_Q["TotalQuality"]["FFQuality"] = FFq[0]
            except (KeyError, TypeError):
                fragment_Q["TotalQuality"]["FFQuality"] = FFq
            json_file.close()

        fragment_Q["system"] = system.readme
        qualities.append(fragment_Q)

    # ---- Sort based on total quality of a simulation
    fragments = ["total", "tails", "headgroup"]

    for sort_based_on in fragments:
        new_qualities = []
        for i in qualities:
            try:
                if i["TotalQuality"][sort_based_on] > 0:
                    new_qualities.append(i)
            except (KeyError, TypeError):
                continue

        sorted_qualities = sorted(new_qualities, key=lambda i: i["TotalQuality"][sort_based_on], reverse=True)

        outputfile = os.path.join(FMDL_DATA_PATH, "Ranking", "SYSTEM_" + sort_based_on + "_Ranking.json")
        with open(outputfile, "w") as fp:
            json.dump(sorted_qualities, fp, default=str, cls=CompactJSONEncoder)
        print(f"Sorted based on {sort_based_on} quality and saved to {outputfile}")

    new_qualities = []
    for i in qualities:
        try:
            if i["TotalQuality"]["FFQuality"] > 0:
                new_qualities.append(i)
        except (KeyError, TypeError):
            continue

    sorted_qualities = sorted(new_qualities, key=lambda i: i["TotalQuality"]["FFQuality"])

    outputfile = os.path.join(FMDL_DATA_PATH, "Ranking", "SYSTEM_FormFactor_Ranking.json")
    with open(outputfile, "w") as fp:
        json.dump(sorted_qualities, fp, default=str, cls=CompactJSONEncoder)
    print("Sorted based on form factor quality and saved to", outputfile)

    # ---- Sorting best simulations for each lipid
    fragments = ["total", "sn-1", "sn-2", "headgroup"]

    for sort_based_on in fragments:
        for lipid in lipids_set:
            new_qualities = []
            for i in qualities:
                try:
                    if i[lipid][sort_based_on] > 0:
                        new_qualities.append(i)
                except (KeyError, TypeError):
                    continue

            sorted_qualities = sorted(new_qualities, key=lambda i: i[lipid][sort_based_on], reverse=True)

            if sorted_qualities:
                outputfile = os.path.join(FMDL_DATA_PATH, "Ranking", lipid + "_" + sort_based_on + "_Ranking.json")
                with open(outputfile, "w") as fp:
                    json.dump(sorted_qualities, fp, default=str, cls=CompactJSONEncoder)
                print(f"Quality of {sort_based_on} of {lipid} sorted and saved to {outputfile}")


if __name__ == "__main__":
    make_ranking()
