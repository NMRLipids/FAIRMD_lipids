#!/usr/bin/env python3
"""
Creates different types of ranking lists inside the Databank.

It ranks simulations based on their quality against experiments. The ranking lists are stored in
``{FMDL_DATA_PATH}/Ranking/`` folder in CSV format.

**Usage:**

.. code-block:: console

    fmdl_make_ranking

No arguments are needed.
"""

import math
import os

import numpy as np
import pandas as pd

from fairmd.lipids import FMDL_DATA_PATH
from fairmd.lipids.api import get_quality, lipids_set
from fairmd.lipids.core import System, initialize_databank

__all__ = ["make_ranking"]


def _make_composition_string(comp: dict) -> str:
    """Return formatted composition for ranking tables."""
    lips = []
    rats = []
    for k, v in comp.items():
        if k in lipids_set:
            lips += [k]
            rats += [sum(v["COUNT"])]
    # sort lipids alphabetically
    lips = np.array(lips)
    rats = np.array(rats)
    ids = np.argsort(lips)
    lips = lips[ids]
    rats = rats[ids]
    div = math.gcd(*rats)
    comp_s = ":".join(map(str, (rats / div).astype(np.int32)))
    return ":".join(lips) + " (" + comp_s + ")"


def _get_hydration_nan(s: System) -> float:
    """Return hydration (allowing nan)."""
    try:
        hy = np.round(s.get_hydration(), 1)
    except ValueError:
        hy = np.nan  # implicit water
    return hy


def make_ranking() -> None:
    """Make ranking CSV tables."""
    ss = initialize_databank()

    # GLOBAL FF and OP rankings
    res_array = []
    for s in ss:
        record = [
            s["ID"],
            s["FF"],
            _make_composition_string(s["COMPOSITION"]),
            s["TEMPERATURE"],
            _get_hydration_nan(s),
            np.round(get_quality(s, experiment="FF"), 4),
            np.round(get_quality(s, experiment="OP"), 4),
            np.round(get_quality(s, part="tails", experiment="OP"), 4),
            np.round(get_quality(s, part="headgroup", experiment="OP"), 4),
        ]
        res_array.append(record)

    df = pd.DataFrame(res_array)
    df.columns = ["ID", "FF", "composition", "T", "Hydration", "Q_FF", "Q_OP", "Q_OP_tails", "Q_OP_heads"]

    # FF ranking
    qff_df = df[~df["Q_FF"].isna()]
    qff_sorted_df = qff_df.sort_values(by="Q_FF", ascending=False)
    fn = os.path.join(FMDL_DATA_PATH, "Ranking", "FF_ranking.csv")
    qff_sorted_df.to_csv(fn, index=False)

    # OP ranking
    qff_df = df[~df["Q_OP"].isna()]
    qff_sorted_df = qff_df.sort_values(by="Q_OP", ascending=False)
    fn = os.path.join(FMDL_DATA_PATH, "Ranking", "OP_ranking.csv")
    qff_sorted_df.to_csv(fn, index=False)

    # Individual lipid NMR
    res_array = {}
    for s in ss:
        for lip in s.lipids:
            qq = get_quality(s, lipid=lip, experiment="OP")
            if np.isnan(qq):
                continue
            try:
                hq = get_quality(s, lipid=lip, part="headgroup", experiment="OP")
            except KeyError:
                hq = np.nan  # Lipid doesn't have head
            try:
                tq = get_quality(s, lipid=lip, part="tails", experiment="OP")
            except KeyError:
                tq = np.nan  # Lipid doesn't have tails

            record = [
                s["ID"],
                s["FF"],
                _make_composition_string(s["COMPOSITION"]),
                s["TEMPERATURE"],
                _get_hydration_nan(s),
                np.round(get_quality(s, experiment="FF"), 4),
                np.round(qq, 4),
                np.round(tq, 4),
                np.round(hq, 4),
            ]
            if lip not in res_array:
                res_array[lip] = []
            res_array[lip].append(record)

    # OP ranking
    for lipid_name, dataframe in res_array.items():
        df = pd.DataFrame(dataframe)
        df.columns = ["ID", "FF", "composition", "T", "Hydration", "Q_FF", "Q_OP", "Q_OP_tails", "Q_OP_heads"]
        q_df = df[~df["Q_OP"].isna()]
        q_sorted_df = q_df.sort_values(by="Q_OP", ascending=False)
        fn = os.path.join(FMDL_DATA_PATH, "Ranking", f"{lipid_name}_ranking.csv")
        q_sorted_df.to_csv(fn, index=False)


if __name__ == "__main__":
    make_ranking()
