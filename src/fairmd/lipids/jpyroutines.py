"""
Databank routines designed for calling from inside jupyter notebooks.
"""

import numpy as np
import pandas as pd
from IPython.display import Markdown, display


def showTable(sorted_qualities, quality):  # noqa: N802 (API)
    """
    Shows a table of simulation qualities against experimental data.

    :param sorted_qualities: list of dictionaries to be shown, available
                            in folder ``Data/Ranking/``
    :param quality: should be either ``TotalQuality`` or universal lipid name.
                    First one shows the system total quality.
                    Latter shows the individual lipid quality.
    """
    rounding = ["headgroup", "sn-1", "sn-2", "total", "tails", "FFQuality"]
    quality_table = []
    pd.set_option("display.max_rows", None)
    for i in sorted_qualities:
        stored_to_table = []
        for k, v in i[quality].items():
            if k in rounding:
                if v and v != float("inf") and not np.isnan(v):
                    i[quality][k] = round(float(v), 2)

        stored_to_table = i[quality]
        try:
            stored_to_table["Forcefield"] = i["system"]["FF"]
        except (KeyError, TypeError):
            display(Markdown("**FAILURE:** no FF defined for the system"))
            display(i["system"]["path"])
            continue

        molecules = ""
        mol_numbers = ""
        for lipid in i["system"]["COMPOSITION"]:
            molecules = molecules + lipid + ":"
            mol_numbers = mol_numbers + str(np.sum(i["system"]["COMPOSITION"][lipid]["COUNT"])) + ":"
        stored_to_table["Molecules"] = molecules[:-1]
        stored_to_table["Number of molecules"] = " (" + mol_numbers[:-1] + ")"
        stored_to_table["Temperature"] = i["system"]["TEMPERATURE"]
        stored_to_table["ID"] = i["system"]["ID"]
        quality_table.append(stored_to_table)
    display(pd.json_normalize(quality_table))
