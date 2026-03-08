#!/usr/bin/env python3

from abc import ABC, abstractmethod

import numpy as np
import pandas as pd
from unidecode import unidecode

from fairmd.lipids.api import get_OP
from fairmd.lipids.core import System, initialize_databank
from fairmd.lipids.experiment import ExperimentCollection, OPExperiment


class OPDataError(Exception):
    """Our specific exception"""


class OPDataStorer(ABC):
    @abstractmethod
    def store_to_hdf5(self, hdf_fname: str):
        """Store the record"""

    @abstractmethod
    def prepare_dataframe(self):
        """Prepare dataframe for storing"""

    @property
    @abstractmethod
    def ass_id(self):
        """Get id of assoc object"""

    def _prepare_df_common(self, mol, opdict, err_extractor):
        """Common code for dataframe conditioning for storage"""
        smiles = mol.metadata["bioschema_properties"]["smiles"]
        print(smiles)
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
        df = pd.DataFrame(
            {
                "id": np.full(len(opdict), -1, dtype=np.int64),  # int64, init NaN
                "val": np.full(len(opdict), np.nan),
                "err": np.full(len(opdict), np.nan),
            }
        )
        cur_row = 0
        for id_, uname in smi2uname.items():
            for opdict_row in OPDataStorer.by_c_uname(opdict, uname):
                if cur_row >= len(df):
                    break
                df.at[cur_row, "id"] = id_
                df.at[cur_row, "val"] = opdict_row[0]
                df.at[cur_row, "err"] = err_extractor(opdict_row)
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
    @property
    def ass_id(self):
        return self._e.exp_id

    def store_to_hdf5(self, hdf_fname: str):
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
        group += (
            unidecode(self._e.exp_id)
            .replace("-", "_")  # prefer variable-like strs
            .replace("/", "_")
            .replace(".", "")
        )
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

    def __init__(self, e: OPExperiment, lname: str):
        self._e = e
        self._lname = lname

    def prepare_dataframe(self):
        mol = self._e.lipids[self._lname]
        opdict = self._e.data[self._lname]
        self._prepare_df_common(mol, opdict, err_extractor=lambda x: 0.02 if len(x) == 1 else x[1])


class SimOPDataStorer(OPDataStorer):
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


H5_EXP_MASTER = "exp-op-dataset.h5"


def main_exps() -> None:
    print("Generating OP datasets from experiments.")
    exps = ExperimentCollection.load_from_data("OPExperiment")
    for exp in exps:
        print(exp)
        for lname in exp.data:
            ods = ExpOPDataStorer(exp, lname)
            try:
                ods.prepare_dataframe()
            except OPDataError as e:
                print("ERROR: ", e)
                continue
            ods.store_to_hdf5(H5_EXP_MASTER)


H5_SIMS_MASTER = "sims-op-dataset.h5"


def main_sims() -> None:
    print("Generating OP dataset from simulations.")
    sims = initialize_databank()
    for sim in sims:
        print(sim)
        opdata = get_OP(sim)
        if opdata is None:
            continue
        for lname in opdata:
            if opdata[lname] is None:
                continue
            ods = SimOPDataStorer(sim, lname)
            try:
                ods.prepare_dataframe()
            except OPDataError as e:
                print("ERROR: ", e)
                continue
            ods.store_to_hdf5(H5_SIMS_MASTER)


if __name__ == "__main__":
    main_exps()
    main_sims()
