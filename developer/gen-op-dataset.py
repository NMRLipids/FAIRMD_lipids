#!/usr/bin/env python3

import argparse
import logging
import re
import sys
from abc import ABC, abstractmethod

import numpy as np
import pandas as pd
from unidecode import unidecode

from fairmd.lipids._base import SampleComposition
from fairmd.lipids.api import get_OP
from fairmd.lipids.core import System, initialize_databank
from fairmd.lipids.experiment import ExperimentCollection, OPExperiment
from fairmd.lipids.molecules import Molecule


class OPDataError(Exception):
    """Our specific exception"""


class OPDataStorer(ABC):
    """Abstract class for OP data storage"""

    @abstractmethod
    def prepare_dataframe(self) -> None:
        """Prepare dataframe for storing"""

    @property
    @abstractmethod
    def ass_id(self) -> str:
        """Get id of assoc object"""

    @property
    @abstractmethod
    def sample(self) -> SampleComposition:
        """Return sample instance"""

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

    def get_mcontent(self) -> dict:
        """Generate membrane content dictionary"""
        _mcontent = self.sample.membrane_composition(basis="molar")
        rdict = {"name": [], "inchikey": [], "fraction": []}
        for lname, frac in _mcontent.items():
            k = lname
            ik = self.sample.lipids[lname].metadata["bioschema_properties"]["inChIKey"]
            rdict["name"] += [k]
            rdict["inchikey"] += [ik]
            rdict["fraction"] += [frac]
        return rdict

    def get_DF_attrs(self) -> tuple[dict, dict]:
        """Genererate attributes for OP and Composition dataframes"""
        inchikey = self.sample.lipids[self._lname].metadata["bioschema_properties"]["inChIKey"]
        smiles = self.sample.lipids[self._lname].metadata["bioschema_properties"]["smiles"]
        opt_attr = {
            "inchikey": inchikey,
            "smiles": smiles,
        }
        hydration = self.sample.get_hydration()
        scontent = self.sample.solution_composition(basis="molar")
        smp_attr = {
            "hydration": hydration,
            "solution": ", ".join([f"{k:<25} {v * 100:>6.1f}%" for k, v in sorted(scontent.items())]),
        }
        return opt_attr, smp_attr

    @abstractmethod
    def get_H5_gname(self) -> str:
        """Generate name for the record in the H5 table"""

    def store_to_hdf5(self, hdf_fname: str) -> None:
        mcontent = self.get_mcontent()
        opt_attr, smp_attr = self.get_DF_attrs()
        # store all vars and df to the HDF5 table
        group = self.get_H5_gname()
        with pd.HDFStore(hdf_fname, "a") as store:
            # DataFrame table
            store.put(f"{group}/op_values", self._df, format="table", data_columns=True)
            opdf = pd.DataFrame(mcontent)
            print(opdf)
            store.put(f"{group}/sample_table", opdf, format="table", data_columns=True)
            # Metadata attributes - I
            op_storer = store.get_storer(f"{group}/op_values")
            for k, v in opt_attr.items():
                op_storer.attrs[k] = v
            # -//- II
            sample_storer = store.get_storer(f"{group}/sample_table")
            for k, v in smp_attr.items():
                sample_storer.attrs[k] = v


class ExpOPDataStorer(OPDataStorer):
    """OP data storer for experiments"""

    DEFAULT_EXP_HDFNAME = "exp-op-dataset.h5"
    """Default filename for the experimental dataset"""

    @property
    def ass_id(self) -> str:
        return self._e.exp_id

    @property
    def sample(self) -> SampleComposition:
        return self._e

    def get_H5_gname(self) -> str:
        group = "E"
        group += re.sub(r"[^A-Za-z0-9_]", "_", unidecode(self._e.exp_id))
        group += "__" + self._lname
        return group

    def get_DF_attrs(self) -> tuple[dict, dict]:
        opt_attr, smp_attr = super().get_DF_attrs()
        opt_attr ["fmdl_expid"] = self._e.exp_id
        nmr_method = self._e.metadata.get("NMR", {}).get("METHOD", False)
        if nmr_method:
            opt_attr["nmr_method"] = nmr_method
        smp_attr["temperature"] = self._e["TEMPERATURE"]
        return opt_attr, smp_attr

    def __init__(self, e: OPExperiment, lname: str) -> None:
        """Initialize with experiment object and lipid name"""
        self._e: OPExperiment = e
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

    @property
    def sample(self) -> SampleComposition:
        return self._s

    def prepare_dataframe(self):
        mol = self._s.lipids[self._lname]
        opdict = get_OP(self._s)[self._lname]
        self._prepare_df_common(mol, opdict, err_extractor=lambda x: x[2])

    def get_H5_gname(self) -> str:
        return f"SIM_{self._s['ID']}__{self._lname}"

    def get_mcontent(self):
        retdic = super().get_mcontent()
        retdic["number"] = [0] * len(retdic["name"])
        retdic["asymmetry"] = [0] * len(retdic["name"])
        _simcomp = self._s["COMPOSITION"]
        for i, lname in enumerate(retdic["name"]):
            cnt = _simcomp[lname]["COUNT"]
            if isinstance(cnt, int):
                asm = np.nan
                cnt = [cnt / 2, cnt / 2]
            else:
                asm = cnt[0] / sum(cnt)
            retdic["number"][i] = sum(cnt) / 2
            retdic["asymmetry"][i] = asm
        return retdic

    def get_DF_attrs(self) -> tuple[dict, dict]:
        opt_attr, smp_attr = super().get_DF_attrs()
        opt_attr["fmdl_simid"] = self._s["ID"]
        smp_attr["temperature"] =  self._s["TEMPERATURE"]
        ff_name = self._s.readme.get("FF", False)
        if ff_name:
            smp_attr["ff_name"] = ff_name
        return opt_attr, smp_attr

    def __init__(self, sim: System, lname: str) -> None:
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
