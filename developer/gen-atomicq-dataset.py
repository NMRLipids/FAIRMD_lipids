#!/usr/bin/env python3

import numpy as np
import pandas as pd

from fairmd.lipids.api import get_OP
from fairmd.lipids.core import System, initialize_databank
from fairmd.lipids.experiment import ExperimentCollection, OPExperiment


class OPQDataError(Exception):
    """Our specific exception"""


class OPQDataStorer():


    DEFAULT_OPQ_HDFNAME = "atomic-opq-dataset.h5"

    def __init__(self, s: System, lname: str) -> None:
        self._s = s
        self._lname = lname
        self._exp_opdicts: list[pd.DataFrame] = []

    def _cvt_op_df(self, opdict: dict, err_pos: int) -> pd.DataFrame:
        op_clean = pd.DataFrame(
            columns=["c", "h", "val", "err"],
        )
        smi2uname = self._get_smi2uname()
        for smid, uname in smi2uname.items():
            _c_dict = {
                k: [v[0], OPExperiment.DEFAULT_ERROR if err_pos >= len(v) else v[err_pos]]  # check
                for k, v in opdict.items()
                if k.split()[0] == uname
            }
            _c_dict_len = len(_c_dict)
            if _c_dict_len == 0:
                continue
            vearr = np.array(list(_c_dict.values()))
            if _c_dict_len == 1:  # one H
                op_clean.loc[len(op_clean)] = [smid, 1, vearr[0, 0], vearr[0, 1]]
            elif _c_dict_len == 3:  # three H. They are always symmetric.
                op_clean.loc[len(op_clean)] = [smid, 1, np.mean(vearr[:, 0]), np.mean(vearr[:, 1])]
            elif _c_dict_len == 2:  # two H. They could be asymmetric.
                vearr = vearr[np.argsort(np.abs(vearr[:, 0]))]
                op_clean.loc[len(op_clean)] = [smid, 1, vearr[0, 0], vearr[0, 1]]
                op_clean.loc[len(op_clean)] = [smid, 2, vearr[1, 0], vearr[1, 1]]
            else:
                msg = (
                    f"Unexpected number of H for {uname} in "
                    f"instance {self.ass_id} // {self._lname}: {_c_dict_len}."
                    " Cannot store."
                )
                raise OPQDataError(msg)
        return op_clean.astype({"c": np.int64, "h": np.int64, "val": np.float64, "err": np.float64})

    def prepare_sim_dataframe(self) -> None:
        """Prepare dataframe for storing."""
        opdict = get_OP(self._s)[self._lname]
        self._sim_op = self._cvt_op_df(opdict, 2)

    @property
    def ass_id(self) -> str:
        """Get id of assoc object"""
        return self._s["ID"]

    def _get_smi2uname(self) -> dict:
        mol = self._s.lipids[self._lname]
        s2u = {}
        for uname, aprops in mol.mapping_dict.items():
            if "SMILEIDX" in aprops:
                smid = int(aprops["SMILEIDX"])
                s2u[smid] = uname
        if not s2u:
            # NO SMILEIDX. Cannot store.
            msg = f"Instance {self.ass_id} // {self._lname} cannot be stored: we don't have SMILEIDX."
            raise OPQDataError(msg)
        return dict(sorted(s2u.items()))

    def add_experiment_data(self, exp_opdict: dict) -> None:
        """Add experimental OP data to the storer. We will use it for Q estimation."""
        self._exp_opdicts.append(self._cvt_op_df(exp_opdict, 1))

    def average_experiment_data(self) -> None:
        """Average experimental OP data if we have more than one."""
        concdf = pd.concat(self._exp_opdicts, ignore_index=True)
        self._exp_opdict_one: pd.DataFrame = concdf.groupby(["c", "h"]).mean().reset_index()

    def compute_q_points(self) -> None:
        """Compute mean(exp) datapoints and inherit sign from simulation OP value"""
        df1 = self._sim_op.merge(
            self._exp_opdict_one,
            on=["c", "h"],
            how="inner",
            suffixes=("_s", "_e"),
        )
        df1["val_e"] = df1["val_e"].abs() * np.sign(df1["val_s"])
        self._qpoints = df1

    def store_to_hdf5(self, hdf_fname: str) -> None:
        """Store the record"""
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
            store.put(f"{group}/op_values", self._qpoints, format="table", data_columns=True)
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




def load_sims() -> None:
    print("Generating OP dataset from simulations.")
    exps = ExperimentCollection.load_from_data("OPExperiment")
    sims = initialize_databank()
    for sim in sims:
        paired_opedict = sim["EXPERIMENT"].get("ORDERPARAMETER", {})
        for lname in sim.lipids:
            if len(paired_opedict.get(lname, [])) == 0:
                continue
            print(sim)

            ods = OPQDataStorer(sim, lname)
            try:
                ods.prepare_sim_dataframe()
            except OPQDataError as e:
                print("ERROR: ", e)
                continue

            for expid in paired_opedict[lname]:
                _exp = exps.loc(expid)
                ods.add_experiment_data(_exp.data[lname])
            ods.average_experiment_data()

            ods.compute_q_points()
            ods.store_to_hdf5(OPQDataStorer.DEFAULT_OPQ_HDFNAME)


if __name__ == "__main__":
    load_sims()
