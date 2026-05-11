#!/usr/bin/env python3

import argparse
import logging
import re
import sys
from abc import ABC, abstractmethod

import numpy as np
import pandas as pd
from unidecode import unidecode

from fairmd.lipids.api import get_OP
from fairmd.lipids.core import System, initialize_databank
from fairmd.lipids.experiment import ExperimentCollection, OPExperiment
from fairmd.lipids.molecules import Molecule


class OPDataError(Exception):
    """Our specific exception"""


class OPDataStorer(ABC):
    """Abstract class for OP data storage"""

    @abstractmethod
    def store_to_hdf5(self, hdf_fname: str) -> None:
        """Store the record"""

    @abstractmethod
    def prepare_dataframe(self) -> None:
        """Prepare dataframe for storing"""

    @property
    @abstractmethod
    def ass_id(self) -> str:
        """Get id of assoc object"""

    def _prepare_df_common(self, mol: Molecule, opdict: dict, err_extractor: callable) -> None:
        """Condition the dataframe for storage"""
        smi2uname = {}
        for uname, aprops in mol.mapping_dict.items():
            if "SMILEIDX" in aprops:
                smid = int(aprops["SMILEIDX"])
                smi2uname[smid] = uname
        if not smi2uname:
            # NO SMILEIDX. Cannot store.
            msg = f"Instance {self.ass_id} // {self._lname} cannot be stored: we don't have SMILEIDX."
            raise OPDataError(msg)
        smi2uname = dict(sorted(smi2uname.items()))
        df = pd.DataFrame(columns=["id", "val", "err"])
        cur_row = 0
        for id_, uname in smi2uname.items():
            for opdict_row in OPDataStorer.by_c_uname(opdict, uname):
                df.loc[cur_row] = [id_, opdict_row[0], err_extractor(opdict_row)]
                cur_row += 1
        self._df = df

    @staticmethod
    def by_c_uname(opdict: dict, uname: str) -> dict:
        """Return only OP vals of named heavy atoms. Sorted by val."""
        res = []
        for k, v in opdict.items():
            if k.split()[0] == uname:
                res.append(v)
        return sorted(res, key=lambda x: x[0])


class ExpOPDataStorer(OPDataStorer):
    """OP data storer for experiments"""

    DEFAULT_EXP_HDFNAME = "exp-op-dataset.h5"
    """Default filename for the experimental dataset"""

    @property
    def ass_id(self) -> str:
        return self._e.exp_id

    def store_to_hdf5(self, hdf_fname: str) -> None:
        _mcontent = self._e.membrane_composition(basis="molar")
        mcontent = {"name": [], "inchikey": [], "fraction": []}
        for lname, frac in _mcontent.items():
            k = lname
            ik = self._e.lipids[lname].metadata["bioschema_properties"]["inChIKey"]
            mcontent["name"] += [k]
            mcontent["inchikey"] += [ik]
            mcontent["fraction"] += [frac]
        hydration = self._e.get_hydration()
        scontent = self._e.solution_composition(basis="molar")
        temperature = self._e["TEMPERATURE"]
        inchikey = self._e.lipids[self._lname].metadata["bioschema_properties"]["inChIKey"]
        smiles = self._e.lipids[self._lname].metadata["bioschema_properties"]["smiles"]
        nmr_method = self._e.metadata.get("NMR", {}).get("METHOD", False)
        # store all vars and df to the HDF5 table
        group = "E"
        group += re.sub(r"[^A-Za-z0-9_]", "_", unidecode(self._e.exp_id))
        group += "__" + self._lname
        with pd.HDFStore(hdf_fname, "a") as store:
            # DataFrame table
            store.put(f"{group}/op_values", self._df, format="table", data_columns=True)
            store.put(f"{group}/sample_table", pd.DataFrame(mcontent), format="table", data_columns=True)
            # Metadata attributes - I
            op_storer = store.get_storer(f"{group}/op_values")
            upd_attr = {
                "inchikey": inchikey,
                "smiles": smiles,
                "fmdl_expid": self._e.exp_id,
            }
            if nmr_method:
                upd_attr["nmr_method"] = nmr_method
            for k, v in upd_attr.items():
                op_storer.attrs[k] = v
            # -//- II
            sample_storer = store.get_storer(f"{group}/sample_table")
            upd_attr = {
                "temperature": temperature,
                "hydration": hydration,
                "solution": ", ".join([f"{k:<25} {v * 100:>6.1f}%" for k, v in sorted(scontent.items())]),
            }
            for k, v in upd_attr.items():
                sample_storer.attrs[k] = v

    def __init__(self, e: OPExperiment, lname: str) -> None:
        """Initialize with experiment object and lipid name"""
        self._e = e
        self._lname = lname

    def prepare_dataframe(self) -> None:
        """Call dataframe preparation"""
        mol = self._e.lipids[self._lname]
        opdict = self._e.data[self._lname]
        self._prepare_df_common(
            mol,
            opdict,
            err_extractor=lambda x: OPExperiment.DEFAULT_ERROR if len(x) == 1 else x[1],
        )


class SimOPDataStorer(OPDataStorer):
    DEFAULT_SIMS_HDFNAME = "sims-op-dataset.h5"
    """Default Dataset Filename"""

    @property
    def ass_id(self):
        return self._s["ID"]

    def prepare_dataframe(self):
        mol = self._s.lipids[self._lname]
        opdict = get_OP(self._s)[self._lname]
        self._prepare_df_common(mol, opdict, err_extractor=lambda x: x[2])

    def store_to_hdf5(self, hdf_fname: str):
        _mcontent = self._s["COMPOSITION"]
        mcontent = {"name": [], "inchikey": [], "number": [], "asymmetry": []}
        for lname, lip in self._s.lipids.items():
            ik = lip.metadata["bioschema_properties"]["inChIKey"]
            cnt = _mcontent[lname]["COUNT"]
            if isinstance(cnt, int):
                asm = np.nan
                cnt = [cnt / 2, cnt / 2]
            else:
                asm = cnt[0] / sum(cnt)
            mcontent["name"] += [lname]
            mcontent["inchikey"] += [ik]
            mcontent["number"] += [sum(cnt) / 2]
            mcontent["asymmetry"] += [asm]
        hydration = self._s.get_hydration()
        scontent = self._s.solution_composition(basis="molar")
        temperature = self._s["TEMPERATURE"]
        inchikey = self._s.lipids[self._lname].metadata["bioschema_properties"]["inChIKey"]
        smiles = self._s.lipids[self._lname].metadata["bioschema_properties"]["smiles"]
        ff_name = self._s.readme.get("FF", False)
        # store all vars and df to the HDF5 table
        group = f"SIM_{self._s['ID']}__{self._lname}"
        with pd.HDFStore(hdf_fname, "a") as store:
            # DataFrame table
            store.put(f"{group}/op_values", self._df, format="table", data_columns=True)
            store.put(f"{group}/simulation_table", pd.DataFrame(mcontent), format="table", data_columns=True)
            # Metadata attributes - I
            op_storer = store.get_storer(f"{group}/op_values")
            upd_attr = {
                "inchikey": inchikey,
                "smiles": smiles,
                "fmdl_simid": self._s["ID"],
            }
            for k, v in upd_attr.items():
                op_storer.attrs[k] = v
            # -//- II
            sample_storer = store.get_storer(f"{group}/simulation_table")
            upd_attr = {
                "temperature": temperature,
                "hydration": hydration,
                "solution": ", ".join([f"{k:<25} {v * 100:>6.1f}%" for k, v in sorted(scontent.items())]),
            }
            if ff_name:
                upd_attr["ff_name"] = ff_name
            for k, v in upd_attr.items():
                sample_storer.attrs[k] = v

    def __init__(self, sim: System, lname: str):
        self._s = sim
        self._lname = lname


def main_exps(log: logging.Logger) -> tuple[int, int]:
    """Generate dataset for experiments"""
    stat_ok, stat_fail = 0, 0
    log.info("\n\nGenerating OP datasets from experiments.")
    exps = ExperimentCollection.load_from_data("OPExperiment")
    for exp in exps:
        for lname in exp.data:
            log.info("%s // %s", str(exp), lname)
            ods = ExpOPDataStorer(exp, lname)
            try:
                ods.prepare_dataframe()
            except OPDataError as e:
                log.error("[from .prepare_dataframe] %s", str(e))  # noqa: TRY400
                stat_fail += 1
                continue
            else:
                stat_ok += 1
            ods.store_to_hdf5(ExpOPDataStorer.DEFAULT_EXP_HDFNAME)
            log.info("..stored!")
    return stat_ok, stat_fail


def main_sims(log: logging.Logger) -> tuple[int, int]:
    """Generate dataset from simulations"""
    stat_ok, stat_fail = 0, 0
    log.info("\n\nGenerating OP dataset from simulations.")
    sims = initialize_databank()
    for sim in sims:
        opdata = get_OP(sim)
        for lname in opdata:
            if opdata is None or opdata[lname] is None:
                stat_fail += 1
                continue
            log.info("%s // %s", str(sim), lname)
            ods = SimOPDataStorer(sim, lname)
            try:
                ods.prepare_dataframe()
            except OPDataError as e:
                log.error("[from .prepare_dataframe] %s", str(e))  # noqa: TRY400
                stat_fail += 1
                continue
            else:
                stat_ok += 1
            ods.store_to_hdf5(SimOPDataStorer.DEFAULT_SIMS_HDFNAME)
            log.info("..stored!")
    return stat_ok, stat_fail


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="OP Dataset Generator",
        description="""
CI helper script for delivering dataset for Kaggle in HDF5 format.
Two DataFrames are stored for each Sim/Exp-lipid pair:
1. SMILES-aligned values for each atom of the molecule
2. Composition table of membrane part of the system""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--exps", action="store_true", help="Generate DS from experiments")
    parser.add_argument("--sims", action="store_true", help="Generate DS from simulations")
    args = parser.parse_args()

    if not args.exps and not args.sims:
        parser.print_usage()
        sys.exit(1)

    lg = logging.getLogger("cli")
    lg.setLevel(logging.INFO)
    h_stderr = logging.StreamHandler(sys.stderr)
    h_stderr.setLevel(logging.INFO)
    lg.addHandler(h_stderr)

    if args.exps:
        e_ok, e_fail = main_exps(lg)
    if args.sims:
        s_ok, s_fail = main_sims(lg)

    lg.info("=======   STATISTICS   ========")
    if args.exps:
        lg.info(f"Stored experiment-lipid pairs: {e_ok} // failed: {e_fail}")
    if args.sims:
        lg.info(f"Stored simulation-lipid pairs: {s_ok} // failed: {s_fail}")
