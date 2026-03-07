#!/usr/bin/env python3

import numpy as np
import pandas as pd
from unidecode import unidecode

from fairmd.lipids.experiment import ExperimentCollection, OPExperiment

H5_MASTER = "exp-op-dataset.h5"


class OPDataError(Exception):
    """Our specific exception"""


class OPDataStorer:
    def store_to_hdf5(self):
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
        with pd.HDFStore(H5_MASTER, "a") as store:
            # DataFrame table
            store.put(f"{group}/op_values", self._df, format="table", data_columns=True)
            store.put(f"{group}/sample_table", pd.DataFrame(mcontent), format="table", data_columns=True)
            # Metadata attributes
            upd_attr = {
                "inchikey": inchikey,
                "smiles": smiles,
                "fmdl_expid": self._e.exp_id,
            }
            if nmr_method:
                upd_attr["nmr_method"] = nmr_method
            op_storer = store.get_storer(f"{group}/op_values")
            for k, v in upd_attr.items():
                op_storer.attrs[k] = v

            # Metadata attributes
            upd_attr = {
                "temperature": temperature,
                "hydration": hydration,
                "solution": ", ".join([f"{k:<25} {v * 100:>6.1f}%" for k, v in sorted(scontent.items())]),
            }
            sample_storer = store.get_storer(f"{group}/sample_table")
            for k, v in upd_attr.items():
                sample_storer.attrs[k] = v

    @staticmethod
    def by_c_uname(opdict: dict, uname: str) -> dict:
        """Return only OP vals of named heavy atoms. Sorted by val."""
        res = []
        for k, v in opdict.items():
            if k.split()[0] == uname:
                res.append(v)
        return sorted(res, key=lambda x: x[0])

    def __init__(self, e: OPExperiment, lname: str):
        self._e = e
        self._lname = lname

    def prepare_dataframe(self):
        mol = self._e.lipids[self._lname]
        opdict = self._e.data[self._lname]
        smiles = mol.metadata["bioschema_properties"]["smiles"]
        print(smiles)
        smi2uname = {}
        for uname, aprops in mol.mapping_dict.items():
            if "SMILEIDX" in aprops:
                smid = int(aprops["SMILEIDX"])
                smi2uname[smid] = uname
        if not smi2uname:
            # NO SMILEIDX. Cannot store.
            msg = f"Experiment {self._e.exp_id} // {self._lname} cannot be stored: we don't have SMILEIDX."
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
            for opdict_row in self.by_c_uname(opdict, uname):
                if cur_row >= len(df):
                    break
                df.at[cur_row, "id"] = id_
                df.at[cur_row, "val"] = opdict_row[0]
                df.at[cur_row, "err"] = 0.02 if len(opdict_row) == 1 else opdict_row[1]
                cur_row += 1
        self._df = df


def main() -> None:
    print("Generating OP datasets.")
    exps = ExperimentCollection.load_from_data("OPExperiment")
    for exp in exps:
        print(exp)
        for lname in exp.data:
            ods = OPDataStorer(exp, lname)
            try:
                ods.prepare_dataframe()
            except OPDataError as e:
                print("ERROR: ", e)
                continue
            ods.store_to_hdf5()


if __name__ == "__main__":
    main()
